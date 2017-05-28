%test
clc, clear, close all
Data_Loader

burst_hw_1 = burst_hw_matrix(:,1) ; % Load the data of the first one into a var
burst_hw_1_msub = burst_hw_1 - mean(burst_hw_1)  ; % Mean subtract the signal

u_mean = mean(burst_hw_1_msub.^2)