% William Page (587000) - Kevin Rassool (540733)   ;
% Semester 2 2017 - University of Melbourne        ; Started:     15/5/17
% MCEN90018 - Advanced Fluid Dynamics              ; Last Edited: 28/5/17
% Assignment 3
%
% Main Script
%% Import Data
clc, clear, close all

fid = fopen('MATLAB/Data/u_hf_ypos3.bin', 'r');
hf_Y3 = fread(fid, '*float') ;

fid_y = fopen('MATLAB/Data/y.txt','r');
data_y = fscanf(fid_y, '%f');

%% Things that are mostly constant

Re_tau = 14000 ; % Reynolds shear stress
n_pos = 40     ; % Number of wall normal position
Fs = 10e3      ; % Sampling frequency

delta = 0.326  ; % Boundary layer thickness (m)

N = length(hf_Y3) ; % Length of the clipped time series signal
tf = 30        ; % Experiment time (s)
dt = 1/Fs      ; % Time interval
df = 1/(N.*dt) ; % Frequency interval
n_spat  = 0:1:(N/2) ; % All mode numbers up to nyquist
f  = n_spat.*df     ; % Frequency vector to match G/A

%% Low pass / High pass filtering

close all ; % Clear any existing figures

% Implement fourier transformation
Glpf = fft(hf_Y3)./N           ; % Take an FFT of the data, normalise to length
A = sqrt(4*(Glpf.*conj(Glpf))) ; % Amplitude function

hf_f_lim = 100        ; % Freqeuncies beyond which the hotfilm is not reliable
cutoff_hf = hf_f_lim  ; % High frequency cut off

[~,n_c_lpf] = min(abs(f-cutoff_hf)) ; % Find frequencies closest to cut-off
Glpf(n_c_lpf:end-n_c_lpf+2) = 0        ; % Cut off them high frequncies
u_lpf = N.*real(ifft(Glpf))            ; % re-construct the fourier signal

figure ; plot(hf_Y3) ; hold on ; plot(u_lpf) ; axis([0,1e3,0,3e-3]) ;
title('High pass filtered data')

%% High pass filter

Ghpf = fft(hf_Y3)./N  ; % Take an FFT of the data, normalise to length
cutoff_lf   = 5                        ; % Low frequency cut off
[~,n_c_hpf] = min(abs(f-cutoff_lf))   ; % Find frequencies closest to cut-off
Ghpf(1:n_c_hpf)         = 0            ; % Cut off them pre nyqist low frequncies
Ghpf(end-n_c_hpf+2:end) = 0            ; % Cut off them post nyquist low frequncies

u_lpf_hpf = N.*real(ifft(Ghpf))        ; % Re-construct the fourier signal

figure ; plot(hf_Y3) ; hold on ; plot(u_lpf_hpf) ; axis([0,1e3,0,3e-3]) ;
title('Original and high pass filtered data') ; xlabel('Time axis'); ylabel('');

%% Cross correlation manually



%%  Plot some information about the signals
% figure ; plot(up2nyq(2:end),A(up2nyq(2:end))) ; title('Energy Information');
% ylabel('Fourier Amplitude') ; xlabel('Fourier Mode')
%
% figure ; plot(f(2:end),A(up2nyq(2:end))) ; title('Energy Information');
% ylabel('Fourier Amplitude') ; xlabel('Frequency')
%
% axis([0,3*cutof_f,0,6e-5]); % look at up2Hz frequencies
%
% f_over = find (f > cutof_f) ; % Indicies where the frequency is unreliable
% A_sensible = A(1:1:(N/2)+1) ; % Frequencies with real information

% f(f_over) = 0 ;              % Remove the entires where the data wasnt good
% f_allover = [f_over, max(f_over)+1:1:max(f_over)+1+length(f_over)  ];
% G_lpf = Glpf() ;
% A_sensible(f_over) = 0 ;     % Remove the amplitudes where the data wasnt good
% A(f_over) = 0 ; % Remove the amplitudes where the data wasnt good

% fr = [f,flip(f)]; % Reconstructed frequency signal
% Ar = [A_sensible;flip(A_sensible)].'; % reconstructed amplitude signal

% u_lpf = real(ifft(G_lpf)) ; % re-construct the signal from the clipped data
% figure ; plot(u_lpf(1:N/2+1)) % Plot the result


%% Qn3, Cross Correlation
%%%%MUST USE FILTERED DATA!!!! %%%

fid_hw = fopen('MATLAB/Data/u_hw_ypos3.bin', 'r');
hw_Y3  = fread(fid_hw, '*float') ;
clip=50;

%Choose vectors to compare %%WHICH ONE SHOULD BE HW and WHICH SHOULD BE HF?
search_r = hw_Y3(1:clip) ;  % Template vector
template = hf_Y3(1:clip) ;  % Search region vector

%Mean subtract
template_ms=template-mean(template);
search_r_ms=search_r-mean(search_r);
%std deviation normalisation
template_ms_std=template_ms/std(template_ms);
search_r_ms_std=search_r_ms/std(search_r_ms);

%pad vectors
temp_zp3 = [zeros(size(template_ms));template_ms_std;zeros(size(template_ms))];
sear_zp3 = [zeros(size(search_r_ms));template_ms_std;zeros(size(search_r_ms))];
N=length(template);

%std deviation of search region
std_search_r_ms=std(search_r_ms);

%long hand cross corr
for kk=-N:1:N-1;
    R_lh(kk+N+1)= sum(sear_zp3(N+1:2*N) .* (temp_zp3([N+1:2*N]+kk)));
end
%normalise results
R_lh=R_lh/N;

%plot results
figure;
n_spat=[-N:1:N-1];
plot(n_spat,R_lh)
title('sptail')

%% Compute FFT correlation
temp_zp2 = [zeros(size(template_ms));template_ms_std];
sear_zp2 = [zeros(size(search_r_ms));template_ms_std];
N=length(template);

%calculate cross corr
cc_fft= (1/N)*fftshift(ifft(conj(fft(sear_zp2)).*fft(temp_zp2))); %
%n vector
n_fft=[-N:1:N-1];

figure;
plot(n_fft,cc_fft)
length(cc_fft)
length(n_fft)

%% Convert to required x-axis

Delta_x_fft = -n_fft.*dt.*mean(search_r)   ;   %Mean of hotwire
dx_on_delta=Delta_x_fft/delta;
% Delta_x_spatial = -lag_spatial.*dt.*mean(search_r)   ;   %

%% Plot Cross Corrs

figure;
hold on
plot(search_r,'r-*')
% plot(template,'p-')
legend('sr','template')

% figure;
% plot(Delta_x_fft,cc_fft)
% title('FFT')
% ax = gca; % current axes
% ax.FontSize = 12;
% ax.TickDir = 'out';
% ax.TickLength = [0.02 0.02];
% ax.XLim = [min(Delta_x_spatial) max(Delta_x_spatial)];
% 
% figure;
% plot(Delta_x_spatial,R_st)
% title('Spatial')
% ax = gca; % current axes
% ax.FontSize = 12;
% ax.TickDir = 'out';
% ax.TickLength = [0.02 0.02];
% ax.XLim = [min(Delta_x_spatial) max(Delta_x_spatial)];



%% Compute FFT correlation at every point.
%Load all data in
Data_Loader;
data_size=size(hw_matrix);
N=data_size(1);
n_fft_loop=[-N:1:N-1];

fft_corr_matrix=zeros(data_size(1)*2,data_size(2));
for jj=1:data_size(2)
template=hf_matrix(:,jj);
search_r=hw_matrix(:,jj);
    %Mean subtract
template_ms=template-mean(template);
search_r_ms=search_r-mean(search_r);
%std deviation normalisation
template_ms_std=template_ms/std(template_ms);
search_r_ms_std=search_r_ms/std(search_r_ms);
%zero pad
temp_zp2 = [zeros(size(template_ms));template_ms_std];
sear_zp2 = [zeros(size(search_r_ms));template_ms_std];
    
%calculate cross corr
cc_fft= (1/N)*fftshift(ifft(conj(fft(sear_zp2)).*fft(temp_zp2)));
%determine actual lag variable
fft_corr_matrix(:,jj)=cc_fft';
end
fft_corr_matrix=fft_corr_matrix';

%% Convert to required x-axis



% Delta_x_fft = -lag_fft.*dt.*mean(search_r)   ;   %

