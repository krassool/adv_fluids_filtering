% William Page (587000) - Kevin Rassool (540733)   ;
% Semester 2 2017 - University of Melbourne        ; Started:     15/5/17
% MCEN90018 - Advanced Fluid Dynamics              ; Last Edited: 28/5/17
% Assignment 3 - Spectral Power Density

clc, clear, close all
Data_Loader

% Re-define constants for higher sample rate
Fs     = 30e3       ; % Sampling frequency
tf     = 30         ; % Experiment time (s)
dt     = 1/Fs       ; % Time interval
N      = Fs*tf      ; % Length of the time serires
df     = 1/(N.*dt)  ; % Frequency interval
n_spat = 0:(N/2)    ; % All mode numbers up to nyquist
f      = n_spat.*df ; % Frequency vector to match G/A

wire_burst      = burst_hw_matrix(:,1) ; % Load the data of the first one into a var
burst_hw_1_msub = wire_burst - mean(wire_burst)  ; % Mean subtract the signal

fft_wireb  = fft(burst_hw_1_msub)./N   ; % Take an FFT of the data, normalise to length
fft_wireb  = fft_wireb(1:N/2+1)        ; % Clip off the post nyquist bs
phi        = 2.*fft_wireb.*conj(fft_wireb)./df ; % Define power spectral density
phi        = phi.';                      % Transpose phi for pre-multiplication

A_spectral = trapz(f,phi);              % Area under the power spectra plot
variance   = mean(burst_hw_1_msub.^2);  % Signal Variance

zero_if_same = A_spectral - variance   % Check if signals same
