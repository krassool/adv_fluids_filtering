
A1  = 5*sin(th) + 3 ;
A11 = 5*sin(th) + 6 ;

A2 = 1*sin(2*th);
A3 = 2.5*cos(7*th);
A4 = 2*cos(3*th);

sumS = A1+A2+A3+A4;

GA1  = fft(A1)./N                 ; % Take an FFT of the data, normalise to length
GA11 = fft(A11)./N                ; % Take an FFT of the data, normalise to length

ampA1  = sqrt(4*(GA1.*conj(GA1)))   ;
ampA11 = sqrt(4*(GA11.*conj(GA11))) ;

figure ; plot(ampA1,'kp') ; hold on  ; plot(ampA11,'rp')