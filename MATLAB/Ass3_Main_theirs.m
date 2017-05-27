% William Page (587000) - Kevin Rassool (540733)   ;
% Semester 2 2017 - University of Melbourne        ; Started:     15/5/17
% MCEN90018 - Advanced Fluid Dynamics              ; Last Edited: 28/5/17
% Assignment 3
%
% Main Script
%% Import Data
clc, clear, close all

fid = fopen('MATLAB/Data/u_hf_ypos3.bin', 'r');
hf_Y3 = fread(fid, '*double') ;

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
n  = 0:1:(N/2) ; % All mode numbers up to nyquist
f  = n.*df     ; % Frequency vector to match G/A

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

close all

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
hw_Y3  = fread(fid_hw, '*double') ;

%Choose vectors to compare
template = hw_Y3(1:22) ;  % Template vector
sr       = hw_Y3(1:22) ;  % Search region vector

% Pad the search region with zeros
sr_padded = pad_vector(template,sr); % Create padded search region

% Calculate the std dev of template
std_d_tem = std(template);

%determine the maximum lag and initialise results vector
lags  = length(sr_padded)-length(template) ; 
R_st  = zeros(length(lags)) ;

tem_sym   = template-mean(template)             ; % T tranformed to be symmetric 

tic
for j = 1:lags+1;
    sr_shifted = circshift(sr_padded,-j+1);        % Shift the matrix by the lag
    sr_clip    = sr_shifted(1:length(template));   % Clip the sr to template size
%     A_sym = A_clip - A_clip_mean;                % A tranformed to be symmetric

    sr_clip_mean = mean2(sr_clip);                % Mean of this region
    sr_sym = sr_clip - sr_clip_mean;              % A tranformed to be symmetric 
%               
%         Rnum  = sum(sum((tem_sym).*(sr_sym))); % First factor for R
%         Rden1 = sum(sum((tem_sym).^2));        % Denominator f1
%         Rden2 = sum(sum((sr_sym).^2));         % Denominator f2
%         R     = (Rnum)/(sqrt(Rden1*Rden2));    % Cross correlation coefficient
%         R_st(j) = R;                           % Store correlation coeff

    Rnum  = sum(sr_clip.*(template));      % First factor for R
    R     = (Rnum);%/(std(sr_clip)*std_d_tem);      % Cross correlation coefficient
    R_st(j) = R;                           % Store correlation coeff
end

%determine actual lag variable
lag_spatial=(1:length(R_st))-floor(length(R_st)/2)-1;

time_xcorr_oldschool = toc

%% Compute FFT correlation
%zero pad vectors
temp_zp= pad_vector(template,sr)
sr_zp= pad_vector(sr,template)

%calculate cross corr
cc_fft= fftshift(ifft(conj(fft(temp_zp)).*fft(sr_zp)));
%determine actual lag variable
lag_fft=(1:length(cc_fft))-floor(length(cc_fft)/2)-1;

%% Plot Cross Corrs

figure;
plot(lag_fft,cc_fft)
title('FFT')
ax = gca; % current axes
ax.FontSize = 12;
ax.TickDir = 'out';
ax.TickLength = [0.02 0.02];
ax.XLim = [min(lag_spatial) max(lag_spatial)];

figure;
plot(lag_spatial,R_st)
title('Spatial')
ax = gca; % current axes
ax.FontSize = 12;
ax.TickDir = 'out';
ax.TickLength = [0.02 0.02];
ax.XLim = [min(lag_spatial) max(lag_spatial)];




