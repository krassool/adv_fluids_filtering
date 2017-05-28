% William Page (587000) - Kevin Rassool (540733)   ;
% Semester 2 2017 - University of Melbourne        ; Started:     15/5/17
% MCEN90018 - Advanced Fluid Dynamics              ; Last Edited: 28/5/17
% Assignment 3 - Conditional averaging

close all % Close any open figures

figure ; plot(film20_lpf) ; title('Filtered input signal') ; % Plot input

condition     = 0   ; % Condition on which the signal is averaged
gradient_fn   = film20_lpf(2:end) - film20_lpf(1:end-1) ; % Discrete Gradient
cond_indicies = find(gradient_fn>condition) ; % Indicies where conditions met

virtual_pad   = 480 ; % Length to look in either direction
cond_average  = 0   ; % Condition on which the signal is averaged

% Remove indicies less than 0 of more than the length of the array
cond_indicies( cond_indicies < virtual_pad ) = [] ; 
cond_indicies( cond_indicies > length(cond_indicies)- virtual_pad ) = [] ;

% Add to the ensemble average where the conditions are met
for i = 1:length(cond_indicies)
    idx = cond_indicies(i);
    cond_average = cond_average+wire20_lpf(idx-virtual_pad:idx+virtual_pad);
end
cond_average = cond_average./length(cond_indicies); % Normalise to correct magnitude

% Create a plot
xplotx = linspace(-virtual_pad,virtual_pad,2*virtual_pad+1);
figure ; plot(xplotx,cond_average) ; xlabel('\Delta t (s)'); ylabel('Average velocity (m/s)');
title('Hotwire signal, conditionally averaged for positive hotfilm gradients')  ;
figure_format(1)