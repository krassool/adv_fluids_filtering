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
Re_theta_rough= zeros(size(dplus_E));
Cf_rough       = zeros(size(dplus_E));

%Rough variables
U_tau_rough    = zeros(size(dplus_E));
ks=325*1e-6; %[m]

%fsolve set up
fsolve_opt=optimset('Display', 'off')

for jj = 1:length(dplus_E)
    
    %this function solves for U_tau for each unique value of delta
        Function = @(U_tau) (1/kap)*log(dplus_E(jj)) - (1/(3*kap))+(2*PI/kap) ...
        - (1/kap)*log(ks*U_tau/nu) + 8.5 - U_sub/U_tau;
    
    
    %input first guess value in to function
    init_val = 0.1;
    
    %Take the next inital value as the result from previous Utau function
    if jj ~= 1
        init_val = U_tau_rough(jj-1);
    end
    % Call fucntion and solve for U tau
    U_tau_rough(jj) = fsolve(Function, init_val, fsolve_opt);
    
    
    %define z_plus
    z_plus = logspace(-4,log10(dplus_E(jj)),n_point);
    
    eta = z_plus./dplus_E(jj);

    %calculate ks plus from U_tau value
    ks_plus= ks*U_tau_rough(jj)/nu;
    %U plus
    u_plus = (1/kap) .*log(z_plus) - (1/kap) .*log(ks_plus) + 8.5 ...
        -(1/(3*kap)) .*eta.^3 + (PI/kap) .* 2 .* (eta.^2) .* (3-2.*eta);
    %S 
    S = (1/kap).*log(dplus_E(jj))-(1/kap).*log(ks_plus)+8.5-(1/(3*kap))+(2*PI/kap);
        % Numerically integrate to get Re_theta
    integrand = (u_plus - (u_plus.^2)./S);
    Re_theta_rough(jj) = trapz(z_plus,integrand);
    % Calculate Cf
    Cf_rough(jj) = 2./(S.^2);
end

figure;
loglog(dplus_E,Re_theta_rough)
title('\delta^+ vs Re_{\theta} for a Rough Wall')
ylabel('Re_{\theta}');
xlabel('\delta^+');
figure_format(1) ;


figure;
semilogx(dplus_E,Cf_rough)
title('\delta^+ vs C_f for a Rough Wall')
ylabel('C_f');
xlabel('\delta^+');
figure_format(1) ;



%% Qn 3, Rex calculation @ x=115m

% Calculate Re_x at every point 
Re_x = cumtrapz(Re_theta_rough, 2./Cf_rough);

% Determine C_f and delta for x = 115 m
x = 115                 ;     %m
Re_x_115 = x*U_sub/nu   ;
Cf_115 = interp1q(Re_x', Cf_rough', Re_x_115)  ;
delta_plus_115 = exp(sqrt(2/Cf_115)*kap - A*kap + 1/3 - 2*PI);
delta_115 = delta_plus_115*nu*sqrt(2/Cf_115)/U_sub

figure;
semilogx(Re_x,Cf_rough)
title('Re_x vs C_f for a Smooth Wall')
hold on
ylabel('C_f');
xlabel('Re_x');
plot(Re_x,Cf_115*ones(length(Re_x)),'--','linewidth',2,'color',[255 105 180]./256);
figure_format(1) ;


