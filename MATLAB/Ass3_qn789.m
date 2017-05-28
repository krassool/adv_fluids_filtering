% William Page (587000) - Kevin Rassool (540733)   ;
% Semester 2 2017 - University of Melbourne        ; Started:     15/5/17
% MCEN90018 - Advanced Fluid Dynamics              ; Last Edited: 28/5/17
% Assignment 3
%
% Main Script
%% Import Data
clc, clear, close all

Data_Loader

wire3 = hw_matrix(:,3) ;
film3 = hf_matrix(:,3) ;

wire20 = hw_matrix(:,20) ;
film20 = hf_matrix(:,20) ;

burst_film20 = burst_hf_matrix(:,1) ; 
burst_wire20 = burst_hw_matrix(:,1) ;

%% Things that are mostly constant
    
Re_tau = 14000 ; % Reynolds shear stress
n_pos  = 40    ; % Number of wall normal position
Fs     = 10e3  ; % Sampling frequency
delta = 0.326  ; % Boundary layer thickness (m) 
N = length(hf_matrix) ; % Length of the clipped time series signal
tf = 30               ; % Experiment time (s)
dt = 1/Fs             ; % Time interval
df = 1/(N.*dt)        ; % Frequency interval
n_spat  = 0:1:(N/2)   ; % All mode numbers up to nyquist
f  = n_spat.*df       ; % Frequency vector to match G/A

%% Low pass / High pass filtering

close all ; % Clear any existing figures

% Implement fourier transformation
Glpf = fft(film20)./N           ; % Take an FFT of the data, normalise to length
A = sqrt(4*(Glpf.*conj(Glpf))) ; % Amplitude function

hf_f_lim = 100        ; % Freqeuncies beyond which the hotfilm is not reliable
cutoff_hf = hf_f_lim  ; % High frequency cut off

[~,n_c_lpf] = min(abs(f-cutoff_hf)) ; % Find frequencies closest to cut-off
Glpf(n_c_lpf:end-n_c_lpf+2) = 0     ; % Cut off them high frequncies
wire20_lpf = N.*real(ifft(Glpf))     ; % Re-construct the fourier signal

figure ; plot(film3) ; hold on ; plot(wire20_lpf) ; 
title('Comparison of original and high pass filtered data')

%% Conditional averaging

close all

figure ; plot(wire20_lpf) ; title('Filtered input signal')

condition     = 0   ;
gradient_fn   = wire20_lpf(2:end) - wire20_lpf(1:end-1) ;
cond_indicies = find(gradient_fn>condition) ;

virtual_pad   = 1e2 ;
cond_average  = 0   ;

cond_indicies( cond_indicies < virtual_pad ) = [];
cond_indicies( cond_indicies > length(cond_indicies)- virtual_pad ) = [] ;

for i = 1:length(cond_indicies)
    idx = cond_indicies(i);
    cond_average = cond_average+wire20(idx-virtual_pad:idx+virtual_pad);
end

cond_average = cond_average./length(cond_indicies);
figure ; plot(cond_average) ; 
title('Conditionally averaged plot')

%%  Spectral Power Density
close all

% Things that were mostly constant before but now are different constant numbers

Fs     = 10e3       ; % Sampling frequency
delta  = 0.326      ; % Boundary layer thickness (m) 
tf     = 3          ; % Experiment time (s)
dt     = 1/Fs       ; % Time interval
df     = 1/(N.*dt)  ; % Frequency interval
N      = Fs*tf      ; % Length of the time serires
n_spat = 0:1:(N/2)  ; % All mode numbers up to nyquist
f      = n_spat.*df ; % Frequency vector to match G/A

burst_hw_1 = burst_hw_matrix(:,1) ; % Load the data of the first one into a var
burst_hw_1_msub = burst_hw_1 - mean(burst_hw_1)  ; % Mean subtract the signal

fft_wireb  = fft(burst_hw_1_msub)./N   ; % Take an FFT of the data, normalise to length

A   = sqrt(4*(fft_wireb_msub.*conj(fft_wireb_msub))) ; % Amplitude function

% HEAD CHECK 1 -> PARSEVALS THEORY HOLDS
A2  = 4.*fft_wireb_msub.*conj(fft_wireb_msub)        ;
sA2 = sum(A2)/4
meanu = mean(burst_hw_1_msub.^2)

% HEAD CHECK 2 -> MODIFIED PARSEVALS HOLDS (DOESNT CURRENTLY HOLD!!)
nnn   = sum(A2(1:N/2+1))/2;
A2_o2 = A2(1:((N/2)+1));
sA2_o2 = sum(A2_o2)/2
meanu2 = mean(burst_hw_1_msub.^2)

% CREATE POEWR SPECTRAL DENSITY FUNCTION
A2_half = A2(1:N/2+1).'  ;
phi = (A2_half)./(2*f) ;

figure ; semilogx(f,A2_half);
title('Amplitude Function on Log Axes');

figure ; semilogx(f,phi);
title('Energy spectral density');

A_spec = trapz(f(2:5),phi(2:5))