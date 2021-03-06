% William Page (587000) - Kevin Rassool (540733)   ;
% Semester 2 2017 - University of Melbourne        ; Started:     15/5/17
% MCEN90018 - Advanced Fluid Dynamics              ; Last Edited: 28/5/17
% Assignment 3
%
% Main Script
%% Set-up
clc, clear, close all
format shortE

%% 
knots_2_mps = 463/900 ; % Conversion between knots and m/s
sub_speed   = 30      ; % Submarine speed

syms y Ut nu A DUp kap PIE delta 
syms U Uinf uplus

eta    = (y/delta);

% Define wall normal coordinates
dplus = (delta*Ut)/nu;
yplus = (y*Ut)/nu;

uplus  = (1/kap)*log(yplus) + A - DUp - (eta^3)*(1/(3*kap)) + PIE/(kap) ...
         * (2*(eta^2)) * (3 - 2*eta) ;
     
S      = (1/kap)*log(dplus)+ A - DUp - (1/(3*kap)) + (2*PIE/kap);
dzplus = (nu/Ut)*(uplus/S - (uplus/S)^2 );

to_print = [dplus,yplus,uplus, S , dzplus];
myprint(to_print)

% pretty(U_ofy) % Prints out the equation in an easy to see way

% Define smooth wall constants
A    = 5        ; % Experimental constants
DUp  = 0        ; % Delta U+
kap  = 0.41     ; % Kappa, von karman constant
nu   = 8.97e-7  ; % Kinematic viscoity (m^2 . s^-1)

dplus_E = logspace(1e2,1e6,1e3);



dtheta_smooth    = simplify(subs(dzplus))
U_ofy_smooth_vpa = vpa(dtheta_smooth,4)

pretty(U_ofy_smooth_vpa)
