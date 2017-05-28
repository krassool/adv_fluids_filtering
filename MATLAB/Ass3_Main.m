% William Page (587000) - Kevin Rassool (540733)   ;
% Semester 2 2017 - University of Melbourne        ; Started:     15/5/17
% MCEN90018 - Advanced Fluid Dynamics              ; Last Edited: 28/5/17
% Assignment 3
%
% Main Script
%% Import Data
clc, clear, close all

Data_Loader;
hf_Y3=hf_matrix(:,3);

fid_y = fopen('MATLAB/Data/y.txt','r');
data_y = fscanf(fid_y, '%f')/1000;

%% Things that are mostly constant

Re_tau = 14000 ; % Reynolds shear stress
n_pos = 40     ; % Number of wall normal position
Fs = 10e3      ; % Sampling frequency
delta = 0.326  ; % Boundary layer thickness (m)

N_3 = length(hf_Y3) ; % Length of the clipped time series signal
tf = 30        ; % Experiment time (s)
dt = 1/Fs      ; % Time interval
df = 1/(N_3.*dt) ; % Frequency interval
n_spat  = 0:1:(N_3/2) ; % All mode numbers up to nyquist
f  = n_spat.*df     ; % Frequency vector to match G/A

lpf_cutoff_value = 100 ;
hpf_cutoff_value = 5   ;
%% Low pass Filter
%Call LPF function
u_lpf=Low_Pass_Filter(hf_Y3,lpf_cutoff_value,f);
%plot results
lpffig = figure ; 
hold on ;  plot([1:N_3]*dt,hf_Y3,'g','linewidth',1) ; plot([1:N_3]*dt,u_lpf,'k','linewidth',4); plot([1:N_3]*dt,u_lpf,'linewidth',2,'color',[255 105 180]./256) ; %axis([0,1e3,0,3e-3]) ;
legend('Pre-filter Data','Filtered Data')
title('Low pass filtered data')
xlim([1 2])
xlabel('t [s]')
ylabel('u [m/s]')
figure_format(1) ;

%% High pass filter
%Call HPF function
u_hpf=High_Pass_Filter(hf_Y3,hpf_cutoff_value,f);

%plot results
%high pass filtered with mean subtracted original
lpffig = figure ; 
hold on ;  plot([1:N_3]*dt,hf_Y3-mean(hf_Y3),'g','linewidth',1) ; plot([1:N_3]*dt,u_hpf,'k','linewidth',4); plot([1:N_3]*dt,u_hpf,'linewidth',2,'color',[255 105 180]./256) ; %axis([0,1e3,0,3e-3]) ;
legend('Pre-filter Data','Filtered Data')
title('High pass filtered data')
xlim([1 1.4])
xlabel('t [s]')
ylabel('u [m/s]')
figure_format(1) ;

%% Qn3, Cross Correlation
%%%%MUST USE FILTERED DATA!!!! %%%

clip=N_3; %this should = N_3
hw_Y3=hw_matrix(:,3);

%Choose vectors to compare %%WHICH ONE SHOULD BE HW and WHICH SHOULD BE HF?
search_r = hf_Y3(1:clip) ;  % Serach Region vector is Hot Film
template = hw_Y3(1:clip) ;  % Template vector is Hot Wire

%Filter the data
search_r = Low_Pass_Filter(search_r,lpf_cutoff_value,f); %Filter the data
template = Low_Pass_Filter(template,70,f); %Filter the data

%Filter the data
% search_r = High_Pass_Filter(search_r,hpf_cutoff_value,f); %Filter the data (hpf)
% template = High_Pass_Filter(template,hpf_cutoff_value,f); %Filter the data (lpf)

%Plot the two signals
figure; 
subplot(2,1,1); plot([1:length(template)]*dt,template,'b'); xlim([0 .5]);
title('Filtered Hot Wire and Hot Film Data (Clipped)');
figure_format(1); xlabel('Time [s]'); ylabel('Velocity [m/s]'); legend('Hot Wire Velocity Profile')
subplot(2,1,2); plot([1:length(template)]*dt,search_r,'r'); xlim([0 .5]);
figure_format(1); xlabel('Time [s]'); ylabel('Velocity [m/s]'); legend('Hot Film Velocity Profile')

%Mean subtract
template_ms=template-mean(template);
search_r_ms=search_r-mean(search_r);
%std deviation normalisation
template_ms_std=template_ms/std(template_ms);
search_r_ms_std=search_r_ms/std(search_r_ms);

%pad vectors
temp_zp3 = [zeros(size(template_ms));template_ms_std;zeros(size(template_ms))];
sear_zp3 = [zeros(size(search_r_ms));search_r_ms_std;zeros(size(search_r_ms))];
N_spat=length(template);

%std deviation of search region
std_search_r_ms=std(search_r_ms);

% %%%%%%long hand cross corr %%%%%%%%%%%
% for kk=-N_spat:1:N_spat-1;
%     R_lh(kk+N_spat+1)= sum(sear_zp3(N_spat+1:2*N_spat) .* (temp_zp3([N_spat+1:2*N_spat]+kk)));
% end
% %normalise results
% R_lh=R_lh/N_spat;
% n_spat=[-N_spat:1:N_spat-1];


%% Compute FFT correlation Low Pass Filter
temp_zp2 = [zeros(size(template_ms));template_ms_std];
sear_zp2 = [zeros(size(search_r_ms));search_r_ms_std];
N_fft=length(template);

%calculate cross corr
cc_fft= (1/N_fft)*fftshift(ifft(conj(fft(sear_zp2)).*fft(temp_zp2))); %
%n vector
n_fft=[-N_fft:1:N_fft-1];

%% Convert to required x-axis

Delta_x_fft = -n_fft.*dt.*mean(search_r)   ;   %Mean of hotwire
dx_on_delta_fft =Delta_x_fft/delta;
Delta_x_spatial = -n_spat.*dt.*mean(search_r)   ;   %
dx_on_delta_spat =Delta_x_spatial/delta;

%% Plot Cross Corrs

figure;
hold on
plot(dx_on_delta_fft,cc_fft,'linewidth',2,'color',[255 105 180]./256);    %plot freq
% plot(dx_on_delta_spat,R_lh,'g--','linewidth',2);   %plot spatial
ax = gca; % current axes
ax.FontSize = 12;
ax.TickLength = [0.02 0.02];
ax.XLim = [min(Delta_x_fft) max(Delta_x_fft)];
figure_format(1) ;

% legend()
title('Cross Correlation')
xlabel('\Delta x/ \delta ')
ylabel('R')
legend('Frequency Domain Cross Correlation','Spatial Domain Cross Correlation')

%% Compute FFT correlation at every point.
%Load all data in
data_size=size(hw_matrix);
N=data_size(1);
n_fft_loop=[-N:1:N-1];
%initialise data
fft_corr_matrix=zeros(data_size(1)*2,data_size(2));
dx_on_delta_loop=zeros(length(n_fft_loop),data_size(2));

for jj=1:data_size(2)
template=Low_Pass_Filter( hw_matrix(:,jj),lpf_cutoff_value,f);
search_r=Low_Pass_Filter( hf_matrix(:,jj),lpf_cutoff_value,f);
    %Mean subtract
template_ms=template-mean(template);
search_r_ms=search_r-mean(search_r);
%std deviation normalisation
template_ms_std=template_ms/std(template_ms);
search_r_ms_std=search_r_ms/std(search_r_ms);
%zero pad
temp_zp2 = [zeros(size(template_ms));template_ms_std];
sear_zp2 = [zeros(size(search_r_ms));search_r_ms_std];
    
%calculate cross corr
cc_fft= (1/N)*fftshift(ifft(conj(fft(sear_zp2)).*fft(temp_zp2)));
%determine actual lag variable
fft_corr_matrix(:,jj)=cc_fft';

Delta_x_loop = -n_fft_loop.*dt.*mean(hw_matrix(:,jj))   ;   %Mean of hotwire
dx_loop(:,jj) =Delta_x_loop;
end
fft_corr_matrix=fft_corr_matrix';

%% Plot resulting matrix to on
% 
figure;
pcolor(dx_loop',data_y/delta, fft_corr_matrix); shading interp

set(gca,'FontSize',16)
colorbar
colormap(jet(20))
caxis([0 max(max(fft_corr_matrix))])
xlim([-3 3])
ylim([0 max(data_y/delta)])
xlabel('\Delta x')
ylabel('y/\delta')



%%
% figure;
% for ll=1:data_size(2)
% contour(dx_on_delta_loop',data_y/delta, fft_corr_matrix,10)
% xlim([-3 3])
% ylim([0 max(data_y/delta)])
% CONTOUR WORKS APPROX:
      %contour(dx_on_delta_loop(:,1)',data_y/delta, fft_corr_matrix,10)
%this contour doesn't have the right values for x/delta

% fitfns = figure ; figure_format() ;
% hold on ; plot(u_pr,V_pr) ; plot(u_po,V_po) ;
% legend('Pre-experiment calibration','Post-experiment calibration')
% title('Fitted calibration fns')

