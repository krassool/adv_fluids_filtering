% Assignment 1  - Advanced Fluid Mechanics - MCEN90018 
% William Page SEM1 2017
% QN1 - Calculation of cross correlation coefficient
% Last Edited: 17/3/17

function [tx,ty,R_st] = ass1_q1(t,A) 
%
% Cross correlation coeffient calculation, conditions are:
%       - length(t) <= length(A)
%       - t and A are 2 dimensional arrays 

% [Mt,Nt] = size(t); % Size info for t 
% [MA,NA] = size(A); % Size info for A
% dim_t = ndims(t)  % Dimensions of array t 
% dim_A = ndims(A)  % Dimensions of array A 
% Do some case elimination (later)

SA = size(A) ; St = size(t)       ; % Size of input functions       
A_pad       = zeros(SA+2*St); % maybe change this to 1 later, 0 ng for xcorr
A_pad(St(1)+1:St(1)+SA(1) , St(2)+1:St(2)+SA(2)) = A ;
figure ; spy(A_pad) ;

lags        = size(A_pad)-size(t) ; % Maximum lag values
tmean       = mean2(t)            ; % Mean of input signal t
t_sym       = t-tmean             ; % t tranformed to be symmetric 
R_st        = zeros(lags) ;  % Pre-allocate cross corr array

tic
for i = 1:lags(1)
    for j = 1:lags(2)
        SRx = i:i+(length(t(:,1))-1);            % Search range in A_Look
        SRy = j:j+(length(t(1,:))-1);            % Search range in A_Look

        A_clip = A_pad(SRx,SRy);                 % Clipped A to match lagged t
        A_clip_mean = mean2(A_clip);             % Mean of this region

        A_sym = A_clip - A_clip_mean;            % A tranformed to be symmetric 
              
        Rnum  = sum(sum((t_sym).*(A_sym)));      % First factor for R
        Rden1 = sum(sum((t_sym).^2));            % denominator f1
        Rden2 = sum(sum((A_sym).^2));            % denominator f2
        R     = (Rnum)/(sqrt(Rden1*Rden2));      % Cross correlation coefficient
        R_st(i,j) = R;                           % Store correlation coeff
    end
end

toc

figure; R_plt = surf(R_st) ; title('my norm xcorr ; R vs Lag') ; 
xlabel('Tau X') ; ylabel('Tau Y') ;
set(R_plt, 'edgecolor','none');

R_matlab = normxcorr2(t,A);

size_mnxc2 = size(R_st)
size_nxc2 = size(R_matlab)

difference = R_st-R_matlab

figure; R_matlab_plt = surf(R_matlab) ; title('norm xcorr ; R vs Lag') ; 
xlabel('Tau X') ; ylabel('Tau Y') ;
set(R_matlab_plt, 'edgecolor','none');

bestFit = max(max(abs(R_st)))
[lag_x,lag_y]  = find(abs(R_st)==bestFit)

[tx,ty] = [lag_x,lag_y]