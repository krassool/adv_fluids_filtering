% William Page (587000) - Kevin Rassool (540733)   ;
% Semester 2 2017 - University of Melbourne        ; Started:     15/5/17
% MCEN90018 - Advanced Fluid Dynamics              ; Last Edited: 28/5/17
% Assignment 3
%
% Main Script
%% Import Data
clc, clear, close all

fid = fopen('MATLAB/Data/u_hf_ypos1.bin', 'r');
data = fread(fid, '*double');

%% See whats in the box today
close all ; % Clear any existing figures

up2  = length(data)          ; % Index up to which we will look at spectra 
clip = data(1:up2)  ; % Implement data clipping

G = fft(clip) ; % Take an FFT of the data
N = length(clip) ; % Length of the clipped time series signal
A = sqrt (4*(G./N).*(conj(G/N)) ) ; % Amplitude 

Fs = 10e3 ; % Sampling frequency
dt = 1/Fs ; % Time interval
df = 1/(N.*dt) ; % Frequency interval
n  = 0:1:(N/2) ; % All mode numbers up to nyquist
f  = n.*df ; % Frequency vector to match G/A

cutof_f = 100 ; % Hz at which the data isnt good 

up2nyq = 1:1:N/2+1 ; % Frequency data that is valid
dt     = 1/20      ; % Sampling interval = 1/f = 10kHz 

figure ; plot(up2nyq(2:end),A(up2nyq(2:end))) ; title('Energy Information'); 
ylabel('Fourier Amplitude') ; xlabel('Fourier Mode')

figure ; plot(f(2:end),A(up2nyq(2:end))) ; title('Energy Information'); 
ylabel('Fourier Amplitude') ; xlabel('Frequency')

axis([0,3*cutof_f,0,6e-5]); % look at up2Hz frequencies

f_over = find (f > cutof_f) ; % Indicies where the frequency is unreliable 
A_sensible = A(1:1:(N/2)+1) ; % Frequencies with real information

f(f_over) = 0 ;              % Remove the entires where the data wasnt good
A_sensible(f_over) = 0;      % Remove the amplitudes where the data wasnt good
A(f_over) = 0;      % Remove the amplitudes where the data wasnt good

fr = [f,flip(f)]; % Reconstructed frequency signal 
Ar = [A_sensible;flip(A_sensible)].'; % reconstructed amplitude signal

u_lpf = real(ifft(Ar)) ; % re-construct the signal from the clipped data
figure ; plot(u_lpf(1:N/2+1)) % Plot the result
