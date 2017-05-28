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

% Do the plots

<<<<<<< HEAD
figure ; semilogx(f,f.*phi) ;
    title('Pre-multiplied Power Spectral Density Plot ')
=======
figure ; semilogx(f,f.*phi) ; 
title('Spectral Power Density Plot')
>>>>>>> 9b2d1b5bb72d8f5645b259c7def3eeb86e7a6fb8
xlabel('frequency - f (Hz)')
ylabel('Freqency \times Power Spectral Density - f \phi')
% figure_format(1)
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
<<<<<<< HEAD
    n_spat_clip   = 0:(N_clip/2)    ; % All mode numbers up to nyquist
    f_clip        = n_spat_clip.*df_clip;  % Frequency vector to match clipped 'g'
    
=======
    n_spat_clip   = 0:N_clip/2;
    f_clip        = n_spat_clip.*df_clip ; % Frequency vector to match clipped 'g'
   
>>>>>>> 9b2d1b5bb72d8f5645b259c7def3eeb86e7a6fb8
    % Implement sprectal analysis
    msub_clip    = wire_bst_clip-mean(wire_bst_clip);
    fft_bst_clip = fft(msub_clip)./N_clip   ; % Take an FFT of the data, normalise to length
    fft_bst_clip = fft_bst_clip(1:N_clip/2+1)    ; % Clip off the post nyquist bs
    phi_clip      = 2.*fft_bst_clip.*conj(fft_bst_clip)./df_clip ; % Define power spectral density
    phi_clip      = phi_clip.';                      % Transpose phi for pre-multiplication
    
<<<<<<< HEAD
    % Figure
=======
    size(f_clip)
    size(phi_clip)
    % Figure 
>>>>>>> 9b2d1b5bb72d8f5645b259c7def3eeb86e7a6fb8
    figure; semilogx(f_clip,f_clip.*phi_clip)
    title('Pre-multiplied Power Spectral Density Plot ')
    xlabel('frequency - f (Hz)')
    ylabel('Freqency \times Power Spectral Density - f \phi')
    figure_format(2);
end

%% ENSEBLE AVERAGE POWER SIGNALS

close all

figure ; semilogx(f,f.*phi) ; 

t_selected     = 5 ;
n_sections     = tf/t_selected  ;
hw_bst_clipmat = zeros(length(burst_hw_matrix)/n_sections,n_sections) ;
t1_wire_bst    = burst_hw_matrix(:,1) ;

% For each row, calcuate the PSD and add it to an average
N_sel          = Fs.*t_selected ; % Number of entries in time series
df_sel         = 1/(N_sel.*dt)  ; % Frequency interval
n_spat_sel     = 0:N_sel/2;
f_sel          = n_spat_sel.*df_sel ; % Frequency vector to match clipped 'g'

PSD_ensemble_1sig   = zeros(N_sel/2+1,1)   ; % Preallocate array

% Make into a matrix
for loopvar1 = 1 : n_sections
    st_var  = N*(loopvar1-1)/n_sections+1;
    end_var = loopvar1*N/n_sections;
    hw_bst_clipmat(:,loopvar1) = t1_wire_bst(st_var:end_var);
end

% Loop through all the clipped indicies and average them
for lv2 = 1:size(hw_bst_clipmat,2)
    tmp_wire = hw_bst_clipmat(:,lv2)        ; % Assign the col we're lookin at
    msub_tmp = tmp_wire - mean(tmp_wire)    ; % Mean subtract the signal
    
    fft_bst_tmp = fft(msub_tmp)./N_sel      ; % Take an FFT of the data, normalise to length
    fft_bst_tmp = fft_bst_tmp(1:N_sel/2+1)  ; % Clip off the post nyquist bs
    
    phi_tmp    = 2.*fft_bst_tmp.*conj(fft_bst_tmp)./df_sel ; % Define power spectral density
    
    PSD_ensemble_1sig = PSD_ensemble_1sig+phi_tmp;
end

PSD_ensemble_1sig = PSD_ensemble_1sig/n_sections;
f_sel = f_sel.';

figure ; semilogx(f_sel,f_sel.*PSD_ensemble_1sig)
title('Second stage of averaging/converging')

% ENSEMBLE AVERAGE THEM FOR EVERY SIGNAL
n_sigs = 10;
PSD_ensemble_allsigs = zeros(N_sel/2+1,1);

for sig = 1:n_sigs
    
    tmp_wire_bst = burst_hw_matrix(:,sig) ;
    
    % Pre-allocate temp variable to store 5s clips in the columns
    hw_bst_clipmat_tmp = zeros(length(burst_hw_matrix)/n_sections,n_sections) ;
    
    % Make a matrix, with each clipped bit in the columns
    % COULD MAKE THIS MUCH MORE ELEGANT USING RESHAPE
    
    for lv_clip1 = 1 : n_sections
        st_var  = N*(lv_clip1-1)/n_sections+1;
        end_var = lv_clip1*N/n_sections;
        hw_bst_clipmat_tmp(:,lv_clip1) = tmp_wire_bst(st_var:end_var);
        
    end
    
    % Loop through all the clipped indicies and average them
    for lv_col = 1:size(hw_bst_clipmat_tmp,2)
        tmp_wire = hw_bst_clipmat_tmp(:,lv_col)  ; % Assign the col we're lookin at
        msub_tmp = tmp_wire - mean(tmp_wire)    ; % Mean subtract the signal

        fft_bst_tmp = fft(msub_tmp)./N_sel      ; % Take an FFT of the data, normalise to length
        fft_bst_tmp = fft_bst_tmp(1:N_sel/2+1)  ; % Clip off the post nyquist bs

        phi_tmp    = 2.*fft_bst_tmp.*conj(fft_bst_tmp)./df_sel ; % Define power spectral density

        PSD_ensemble_allsigs = PSD_ensemble_allsigs + phi_tmp ;
    end
end

PSD_ensemble_allsigs = PSD_ensemble_allsigs./(n_sigs*n_sections);
pm_psd = f_sel.*PSD_ensemble_allsigs;
figure ; semilogx(f_sel,pm_psd)
title('Final stage of averaging/converging')

f_sel_CF = log(f_sel);
[xData, yData] = prepareCurveData( f_sel_CF, pm_psd );
ft = fittype( 'poly9' );% Set up fittype and options.
[fitresult, gof] = fit( xData, yData, ft );% Fit model to data.

% Plot fit with data
xrange = linspace(-0.6931,9.2103,1e3);
y_new = fitresult(xrange );
figure ; h = semilogx(exp(xData),yData,'b.') ; hold on ;
semilogx(exp(xrange),y_new)
xlabel('Frequency') ; ylabel('Pre-Multiplied PSD') ; grid on
legend('Spectral Density Scatter','Approximate energy function')
figure_format(1);