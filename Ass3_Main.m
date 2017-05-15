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

up2  = 1e2          ; % Index up to which we will look at spectra 
clip = data(1:up2)  ; % Implement data clipping

G = fft(clip) ; % Take an FFT of the data
N = length(clip) ; % Length of the clipped time series signal
A = sqrt (4*(G./N).*(conj(G/N)) ) ; % Amplitude 

up2nyq = 1:1:N/2+1 ; % Frequency data that is valid
dt     = 1/20    ; % Sampling interval = 1/f = 10kHz 

n = 0:1:(N/2)      ; % All mode numbers up to nyquist
f = n/(N*dt)       ; % Frequency variable corresponding to A

figure ; plot(up2nyq(2:end),A(up2nyq(2:end))) ; title('Energy Data'); 
ylabel('Fourier Amplitude') ; xlabel('Fourier Mode')

figure ; plot(f(2:end),A(up2nyq(2:end))) ; title('Energy Data'); 
ylabel('Fourier Amplitude') ; xlabel('Frequency')

% data(1:100)
% plot(data(1:100))
