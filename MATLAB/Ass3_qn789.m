% William Page (587000) - Kevin Rassool (540733)   ;
% Semester 2 2017 - University of Melbourne        ; Started:     15/5/17
% MCEN90018 - Advanced Fluid Dynamics              ; Last Edited: 28/5/17
% Assignment 3
%
% Main Script
%% Clear Workspace -> Import Data
clc, clear, close all
Data_Loader

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

wire20 = hw_matrix(:,20) ;
film20 = hf_matrix(:,20) ;

% Implement fourier transformation
Glpf = fft(film20)./N           ; % Take an FFT of the data, normalise to length
A = sqrt(4*(Glpf.*conj(Glpf))) ; % Amplitude function

hf_f_lim = 100        ; % Freqeuncies beyond which the hotfilm is not reliable
cutoff_hf = hf_f_lim  ; % High frequency cut off

[~,n_c_lpf] = min(abs(f-cutoff_hf)) ; % Find frequencies closest to cut-off
Glpf(n_c_lpf:end-n_c_lpf+2) = 0     ; % Cut off them high frequncies
wire20_lpf = N.*real(ifft(Glpf))    ; % Re-construct the fourier signal

figure ; plot(film20) ; hold on ; plot(wire20_lpf) ; 
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

% Do the plots

figure ; semilogx(f,f.*phi) ; 
title('Pre-multiplied energy plot (to make meaningful again)')
xlabel('frequency - f (Hz)')
ylabel('Freqency \times Power Spectral Density - f \phi')
figure_format(1)
%% HEAD CHECKS TO CHECK MATHS IS STILL MATHS

% HEAD CHECK 1 -> PARSEVALS THEORY HOLDS
A2  = 4.*fft_wireb.*conj(fft_wireb)        ;
sA2 = sum(A2)/2
meanu = mean(burst_hw_1_msub.^2)

% HEAD CHECK 2 -> MODIFIED PARSEVALS HOLDS 
A2_o2 = A2(1:((N/2)+1));
sA2_o2 = sum(A2_o2)/2
meanu2 = mean(burst_hw_1_msub.^2)

%% Question 9 : De-hairy function 

% Re-define constants as nessesary
Fs     = 30e3       ; % Sampling frequency
tf     = 30         ; % Experiment time (s)
dt     = 1/Fs       ; % Time interval
N      = Fs*tf      ; % Length of the time serires
df     = 1/(N.*dt)  ; % Frequency interval
n_spat = 0:(N/2)    ; % All mode numbers up to nyquist
f      = n_spat.*df ; % Frequency vector to match G/A

wire_bst = burst_hw_matrix(:,1);
%% Find when the informationn dissappears
close all

for t_clip = [1 5 30]
 
    % Re-define the signal and relevant constants
    wire_bst_clip = wire_bst(1:t_clip*Fs) ;
    N_clip        = length(wire_bst_clip) ;       
    df_clip       = 1/(N_clip.*dt)  ; % Frequency interval for clipped signal
    n_spat_clip   = 0:(N_clip/2)    ; % All mode numbers up to nyquist
    f_clip        = n_spat_clip.*df_clip;  % Frequency vector to match clipped 'g'
   
    % Implement sprectal analysis
    msub_clip    = wire_bst_clip-mean(wire_bst_clip);
    fft_bst_clip = fft(msub_clip)./N_clip   ; % Take an FFT of the data, normalise to length
    fft_bst_clip = fft_bst_clip(1:N_clip/2+1)    ; % Clip off the post nyquist bs
    phi_clip      = 2.*fft_bst_clip.*conj(fft_bst_clip)./df_clip ; % Define power spectral density
    phi_clip      = phi_clip.';                      % Transpose phi for pre-multiplication
    
    % Figure 
    figure; semilogx(f_clip,f_clip.*phi_clip)
    figure_format(2);
end

%% CLIP THEM OUT OF EVERY SIG

% ENSEBLE AVERAGE THEM
