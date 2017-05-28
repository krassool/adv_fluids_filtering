% William Page (587000) - Kevin Rassool (540733)   ;
% Semester 2 2017 - University of Melbourne        ; Started:     15/5/17
% MCEN90018 - Advanced Fluid Dynamics              ; Last Edited: 28/5/17
% Assignment 3
%
% Main Script
%% Set-up
clc, clear, close all

%% Speed Conversion
knots_2_mps = 463/900 ; % Conversion between knots and m/s
sub_speed   = 30      ; % Submarine speed
U_sub=sub_speed*knots_2_mps;

% number of z points
n_point=1000;
% Define smooth wall constants
A    = 5        ; % Experimental constants
DUp  = 0        ; % Delta U+
kap  = 0.41     ; % Kappa, von karman constant
nu   = 8.97e-7  ; % Kinematic viscoity (m^2 . s^-1)
PI   = 0.55     ; % For ZPG turbulant boundary layer
%create log-space
dplus_E = logspace(log10(1e2) , log10(1e6), 1e3);
%initilaise vectors
Re_theta = zeros(size(dplus_E));
Cf       = zeros(size(dplus_E));

for i = 1:length(dplus_E)
    
    z_plus = logspace(-4,log10(dplus_E(i)),n_point);
    
    eta = z_plus./dplus_E(i);

    %U plus
    u_plus = (1/kap).*log(z_plus)+A-(1/(3*kap)).*eta.^3 + ...
        (PI/kap).*2.*(eta.^2).*(3-2.*eta);
    %S 
    S = (1/kap).*log(dplus_E(i))+A-(1/(3*kap))+(2*PI/kap);
    % Numerically integrate to get Re_theta
    inetgrand = (u_plus - (u_plus.^2)./S);
    Re_theta(i) = trapz(z_plus,inetgrand);
    % Calculate Cf
    Cf(i) = 2./(S.^2);
end

figure;
loglog(dplus_E,Re_theta)
title('\delta^+ vs Re_{\theta} for a Smooth Wall')
ylabel('Re_{\theta}');
xlabel('\delta^+');
figure_format(1) ;


figure;
semilogx(dplus_E,Cf)
title('\delta^+ vs C_f for a Smooth Wall')
ylabel('C_f');
xlabel('\delta^+');
figure_format(1) ;

%% Qn 3, Rex calculation @ x=115m

% Calculate Re_x at every point 
Re_x = cumtrapz(Re_theta, 2./Cf);

% Determine C_f and delta for x = 115 m
x = 115                 ;     %m
Re_x_115 = x*U_sub/nu   ;
Cf_115 = interp1q(Re_x', Cf', Re_x_115)
delta_plus_115 = exp(sqrt(2/Cf_115)*kap - A*kap + 1/3 - 2*PI);
delta_115 = delta_plus_115*nu*sqrt(2/Cf_115)/U_sub

figure;
semilogx(Re_x,Cf)
title('Re_x vs C_f for a Smooth Wall')
hold on
ylabel('C_f');
xlabel('Re_x');
plot(Re_x,Cf_115*ones(length(Re_x)),'--','linewidth',2,'color',[255 105 180]./256);
figure_format(1) ;


