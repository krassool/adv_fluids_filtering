% MCEN90018:   Advanced Fluid Dynamics - Assignment 3
% ------------------------------------------------------------------------
% Mischka      Kamener    539030                    Last modified: 28/5/16
% Robert       Haberkern  637517
% 
% Script for Part 2

npoints = 1000;
A       = 5;
k       = 0.41;             % kappa
PI      = 0.55;             % For ZPG turbulent boundary layers
nu      = 8.97e-7;          % m^2/s
U_inf   = 30*0.51444444444; % m/s
ks      = 325e-6;           % m

%% Smooth Wall
% Get logarithmic spacing of d_plus
d_plus   = logspace(2,6,npoints);
Re_theta = zeros(size(d_plus));
Cf       = zeros(size(d_plus));

% Numerically evaluate momentum thickness at each d_plus
for i = 1:length(d_plus)
    % Create array to store incremental values of z_plus. Need to start 
    % just above zero to keep value finite.
    z_plus = logspace(-4,log10(d_plus(i)),npoints);
    
    eta = z_plus./d_plus(i);

    % Calculate u_plus (smooth wall)
    u_plus = (1/k).*log(z_plus)+A-(1/(3*k)).*eta.^3 + ...
        (PI/k).*2.*(eta.^2).*(3-2.*eta);
    % Calculate S (smooth wall)
    S = (1/k).*log(d_plus(i))+A-(1/(3*k))+(2*PI/k);
    % Numerically integrate to get Re_theta
    dRe_theta = (u_plus - (u_plus.^2)./S);
    Re_theta(i) = trapz(z_plus,dRe_theta);
    % Calculate Cf
    Cf(i) = 2./(S.^2);
end

% Calculate theta
U_tau = U_inf.*sqrt(Cf./2);
theta = Re_theta.*nu./U_tau;

% Calculate Re_x
Re_x = cumtrapz(Re_theta, 2./Cf);

% Determine C_f and delta for x = 115 m
x = 115; %m
Re_x_115 = x*U_inf/nu;
Cf_115 = interp1q(Re_x', Cf', Re_x_115)
delta_plus_115 = exp(sqrt(2/Cf_115)*k - A*k + 1/3 - 2*PI);
delta_115 = delta_plus_115*nu*sqrt(2/Cf_115)/U_inf

% Plot figures
figure
loglog(d_plus,Re_theta)
ylabel('Re_{\theta}');
xlabel('Re_{\tau} = \delta^+');

figure
semilogx(d_plus, Cf)
ylabel('C_f');
xlabel('Re_{\tau} = \delta^+');

figure
semilogx(Re_x,Cf)
hold on
semilogx([100 Re_x_115], [Cf_115 Cf_115], '--k')
semilogx([Re_x_115 Re_x_115], [0.001 Cf_115], '--k')
hold off
ylabel('C_f');
xlabel('Re_{x}');
legend('Smooth wall');

return

%% Rough Wall
% Get logarithmic spacing of d_plus
Re_theta_r = zeros(size(d_plus));
Cf_r       = zeros(size(d_plus));
U_tau_r    = zeros(size(d_plus));

% Turn off display for fsolve
options = optimset('Display', 'off');

% Calculate Cf for rough wall
for i = 1:length(d_plus)
    % Create function to solve for Ut
    U_tauFun = @(Ut) (1/k)*log(d_plus(i)) - (1/(3*k))+(2*PI/k) ...
        - (1/k)*log(ks*Ut/nu) + 8.5 - U_inf/Ut;
    % Get initial guss for Ut
    guess = 0.1;
    if i ~= 1
        guess = U_tau_r(i-1);
    end
    % Calculate U_tau
    U_tau_r(i) = fsolve(U_tauFun, guess, options);
    
    % Create array to store incremental values of z_plus. Need to start 
    % just above zero to keep value finite.
    z_plus = logspace(-4,log10(d_plus(i)),npoints);
    
    eta = z_plus./d_plus(i);

    % Calculate ks_plus
    ks_plus = ks*U_tau_r(i)/nu;
    % Calculate u_plus (rough wall)
    u_plus = (1/k).*log(z_plus)-(1/k).*log(ks_plus)+8.5 ...
        -(1/(3*k)).*eta.^3 + (PI/k).*2.*(eta.^2).*(3-2.*eta);
    % Calculate S (rough wall)
    S = (1/k).*log(d_plus(i))-(1/k).*log(ks_plus)+8.5 ...
        -(1/(3*k))+(2*PI/k);
    % Numerically integrate to get Re_theta
    dRe_theta = (u_plus - (u_plus.^2)./S);
    Re_theta_r(i) = trapz(z_plus,dRe_theta);
    % Calculate Cf
    Cf_r(i) = 2./(S.^2);
    
end

% Calculate Re_x
Re_x_r = cumtrapz(Re_theta_r, 2./Cf_r);

% Determine C_f and delta for x = 115 m (rough wall)
Cf_115_r = interp1q(Re_x_r', Cf_r', Re_x_115)
delta_plus_115_r = exp(sqrt(2/Cf_115_r)*k + ...
    log(ks*U_inf*sqrt(Cf_115_r/2)/nu) - 8.5*k + 1/3 - 2*PI);
delta_115_r = delta_plus_115_r*nu*sqrt(2/Cf_115_r)/U_inf

% Plot results
figure
semilogx(Re_x, Cf, Re_x_r, Cf_r);
ylabel('C_f');
xlabel('Re_{x}');
legend('Smooth wall', 'k_s = 325\mum');

inPlot = Re_x > 1e4;
inPlot_r = Re_x_r > 1e6;

figure
semilogx(Re_x(inPlot), Cf(inPlot), Re_x_r(inPlot_r), Cf_r(inPlot_r));
hold on
semilogx([10000 Re_x_115], [Cf_115 Cf_115], '--k')
semilogx([Re_x_115 Re_x_115], [0.001 Cf_115], '--k')
semilogx([10000 Re_x_115], [Cf_115_r Cf_115_r], '--r')
semilogx([Re_x_115 Re_x_115], [Cf_115 Cf_115_r], '--r')
hold off
ylabel('C_f');
xlabel('Re_{x}');
legend('Smooth wall', 'k_s = 325\mum');