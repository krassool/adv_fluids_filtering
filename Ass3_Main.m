% William Page (587000) - Kevin Rassool (540733)   ;
% Semester 2 2017 - University of Melbourne        ; Started:     15/5/17
% MCEN90018 - Advanced Fluid Dynamics              ; Last Edited: 28/5/17
% Assignment 3
%
% Main Script
%% Import Data
clc, clear, close all

fid = fopen('MATLAB/Data/u_hf_ypos1.bin', 'r');
data = fread(fid, '*double');

%% See whats in the box today

G = fft(data);
N = length(G);
A = sqrt (4*(G./N).*(conj(G/N)) );

up2nyq = 1:1:N/2+1;

figure ; plot(up2nyq,A(up2nyq)) ; title('Energy>?')

% data(1:100)
% plot(data(1:100))
