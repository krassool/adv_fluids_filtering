function [u_lpf] = Low_Pass_Filter(data,hf_cutoff,f)
%Low PAss Filter
%Takes the input of the data, the cutoff frequency (in hertz) and a
%frequency vector the same length
% Implement fourier transformation
Glpf = fft(data)./length(data)          ; % Take an FFT of the data, normalise to length
A = sqrt(4*(Glpf.*conj(Glpf))) ; % Amplitude function

[~,n_c_lpf] = min(abs(f-hf_cutoff)) ; % Find frequencies closest to cut-off
Glpf(n_c_lpf:end-n_c_lpf+2) = 0        ; % Cut off them high frequncies
u_lpf = length(data).*real(ifft(Glpf))            ; % re-construct the fourier signal
