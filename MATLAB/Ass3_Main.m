% William Page (587000) - Kevin Rassool (540733)   ;
% Semester 2 2017 - University of Melbourne        ; Started:     15/5/17
% MCEN90018 - Advanced Fluid Dynamics              ; Last Edited: 28/5/17
% Assignment 3
%
% Main Script
%% Import Data
clc, clear, close all

fid = fopen('MATLAB/Data/u_hf_ypos3.bin', 'r');
Y3 = fread(fid, '*double') ;

fid_y = fopen('MATLAB/Data/y.txt','r');
data_y = fscanf(fid_y, '%f');

%% Things that are mostly constant

Re_tau = 14000 ; % Reynolds shear stress
n_pos = 40     ; % Number of wall normal position
Fs = 10e3      ; % Sampling frequency

delta = 0.326  ; % Boundary layer thickness (m)

N = length(Y3) ; % Length of the clipped time series signal
tf = 30        ; % Experiment time (s)
dt = 1/Fs      ; % Time interval
df = 1/(N.*dt) ; % Frequency interval
n  = 0:1:(N/2) ; % All mode numbers up to nyquist
f  = n.*df     ; % Frequency vector to match G/A

%% Low pass / High pass filtering

close all ; % Clear any existing figures

% Implement fourier transformation
Glpf = fft(Y3)./N           ; % Take an FFT of the data, normalise to length
A = sqrt(4*(Glpf.*conj(Glpf))) ; % Amplitude function

hf_f_lim = 100        ; % Freqeuncies beyond which the hotfilm is not reliable 
cutoff_hf = hf_f_lim  ; % High frequency cut off

[~,n_c_lpf] = min(abs(f-cutoff_hf)) ; % Find frequencies closest to cut-off 
Glpf(n_c_lpf:end-n_c_lpf+2) = 0        ; % Cut off them high frequncies
u_lpf = N.*real(ifft(Glpf))            ; % re-construct the fourier signal

figure ; plot(Y3) ; hold on ; plot(u_lpf) ; axis([0,1e3,0,3e-3]) ;
title('High pass filtered data')

%% High pass filter 

close all

Ghpf = fft(Y3)./N  ; % Take an FFT of the data, normalise to length
cutoff_lf   = 5                        ; % Low frequency cut off 
[~,n_c_hpf] = min(abs(f-cutoff_lf))   ; % Find frequencies closest to cut-off 
Ghpf(1:n_c_hpf)         = 0            ; % Cut off them pre nyqist low frequncies
Ghpf(end-n_c_hpf+2:end) = 0            ; % Cut off them post nyquist low frequncies

u_lpf_hpf = N.*real(ifft(Ghpf))        ; % Re-construct the fourier signal

figure ; plot(Y3) ; hold on ; plot(u_lpf_hpf) ; axis([0,1e3,0,3e-3]) ;
title('high pass filtered data')

%% Cross correlation manually 

R = 

%%  Plot some information about the signals
figure ; plot(up2nyq(2:end),A(up2nyq(2:end))) ; title('Energy Information'); 
ylabel('Fourier Amplitude') ; xlabel('Fourier Mode')

figure ; plot(f(2:end),A(up2nyq(2:end))) ; title('Energy Information'); 
ylabel('Fourier Amplitude') ; xlabel('Frequency')

axis([0,3*cutof_f,0,6e-5]); % look at up2Hz frequencies

f_over = find (f > cutof_f) ; % Indicies where the frequency is unreliable 
A_sensible = A(1:1:(N/2)+1) ; % Frequencies with real information

f(f_over) = 0 ;              % Remove the entires where the data wasnt good
f_allover = [f_over, max(f_over)+1:1:max(f_over)+1+length(f_over)  ];
G_lpf = Glpf() ;
% A_sensible(f_over) = 0 ;     % Remove the amplitudes where the data wasnt good
% A(f_over) = 0 ; % Remove the amplitudes where the data wasnt good

fr = [f,flip(f)]; % Reconstructed frequency signal 
% Ar = [A_sensible;flip(A_sensible)].'; % reconstructed amplitude signal

u_lpf = real(ifft(G_lpf)) ; % re-construct the signal from the clipped data
figure ; plot(u_lpf(1:N/2+1)) % Plot the result
