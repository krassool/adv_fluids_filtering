% cross crrelation maths
clc, clear, close all

SIGS = [ 0   -0.0833
         0   -0.0833
         0   -0.0833
   -0.0833   -0.0833
    0.9167    0.9167
   -0.0833   -0.0833
   -0.0833   -0.0833
   -0.0833   -0.0833
   -0.0833   -0.0833
   -0.0833   -0.0833
   -0.0833   -0.0833
   -0.0833   -0.0833 ]

% SIGS = [ 0   -0.0833
%          0   -0.0833
%          0   -0.0833
%          0   -0.0833
%          0    0.9167
%    -0.0833   -0.0833
%     0.9167   -0.0833
%    -0.0833   -0.0833
%    -0.0833   -0.0833
%    -0.0833   -0.0833
%    -0.0833   -0.0833
%    -0.0833   -0.0833 ]

t  = SIGS(:,1)
As = SIGS(:,2)

tmean  = mean2(t)
Asmean = mean2(As)

A_cntrd = As - Asmean
t_cntrd = t  - tmean

numerator = sum((A_cntrd).*(t_cntrd))
sig_t = sum(t_cntrd.^2);
sig_A = sum(A_cntrd.^2);
d3 = sqrt(sig_t)*sqrt(sig_A)
% denominator = sqrt( ((As-Asmean).^2) .* ((t-tmean).^2) )

R = numerator./d3