function [u_hpf] = High_Pass_Filter(data,lf_cutoff,f)
%High PAss Filter
%Takes the input of the data, the cutoff frequency (in hertz) and a
%frequency vector the same length
%Fourier Transform
Ghpf = fft(data)./length(data)  ; % Take an FFT of the data, normalise to length
cutoff_lf   = lf_cutoff;                       ; % Low frequency cut off
[~,n_c_hpf] = min(abs(f-cutoff_lf))    ; % Find frequencies closest to cut-off
Ghpf(1:n_c_hpf)         = 0            ; % Cut off them pre nyqist low frequncies
Ghpf(end-n_c_hpf+2:end) = 0            ; % Cut off them post nyquist low frequncies

u_hpf = length(data).*real(ifft(Ghpf))        ; % Re-construct the fourier signal