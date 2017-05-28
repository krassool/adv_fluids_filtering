% MCEN90018:   Advanced Fluid Dynamics - Assignment 3
% ------------------------------------------------------------------------
% Mischka      Kamener    539030                    Last modified: 28/5/16
% Robert       Haberkern  637517
% 
% Script for Part 1

%% Constants
delta = 0.326;   % m
Fs    = 10*10^3; % Hz
T     = 30;      % s
N     = Fs*T;    % Number of samples
y_i   = 0.0002;  % m
y_f   = 0.450;   % m
N_POS = 40;      % Number of data points

%% Question 1

% Read file
f20 = readBin('f',20);

f20_filtered = lpFilter(f20,T,N,50);

% Plot original and filtered signal over 2 seconds
idx = 10000:20000;
figure
hold on
plot(idx/Fs,f20(idx),'LineWidth',1);
plot(idx/Fs,f20_filtered(idx),'LineWidth',2);
hold off
xlabel('t (s)');
ylabel('Velocity (m/s)');
legend('Unfiltered', 'Filtered');
title('Hotfilm signal for Y=20');

%% Question 2

% Read file
f20 = readBin('f',20);
w20 = readBin('w',20);

f20 = lpFilter(f20,T,N,50);

% Find mean
Uc  = mean(w20);

% Adjust for mean and standard deviation
f20 = f20 - mean(f20); f20 = f20./std(f20);
w20 = w20 - mean(w20); w20 = w20./std(w20);

% Zero pad vector to be shifted
w20_p = [zeros(size(w20)); w20; zeros(size(w20))];
R_cc  = zeros(size(-N:N-1))';

tic
% Compute cross correlation. Slide one signal across the other
for n = -N:N-1;
    R_cc(n+N+1) = sum(f20.*w20_p((N+1:2*N)+n)); 
end
toc

% Normalise by number of elements
R_cc = R_cc./N;

% Determine dt values to plot
lim = 3*delta/Uc;
dt = (-N:N-1)./Fs;
dx = -dt*Uc;

% Plot graph
plot(dx(abs(dt)<lim)/delta,R_cc(abs(dt)<lim))
xlabel('\Deltax/\delta');
ylabel('R');

%% Question 3

% Read files
f20 = readBin('f',20);
w20 = readBin('w',20);

f20 = lpFilter(f20,T,N,50);

% Find mean
Uc  = mean(w20);

% Adjust for mean and standard deviation
f20 = f20 - mean(f20); f20 = f20./std(f20);
w20 = w20 - mean(w20); w20 = w20./std(w20);

% Pad with zeros
f20_p = [zeros(size(f20)); f20];
w20_p = [zeros(size(w20)); w20];

tic
% Compute cross correlation using fft (R_fft goes from -N to N-1)
R_fft = fftshift(ifft(conj(fft(f20_p)).*(fft(w20_p))));
toc

% Normalise by number of elements
R_fft = R_fft./N;

% Determine dt values to plot
lim = 3*delta/Uc;
dt = (-N:N-1)./Fs;
dx = -dt*Uc;

% Plot graph
plot(dx(abs(dt)<lim)/delta,R_fft(abs(dt)<lim),'-', ...
    dx(abs(dt)<lim)/delta,R_cc(abs(dt)<lim),'--')
xlabel('\Deltax/\delta');
ylabel('R');
legend('FFT method', 'Long-hand method');

%% Question 4

% Create meshgrid and initialise matrices
[dt, y] = meshgrid((-N:N-1)./Fs, logspace(log10(y_i), log10(y_f), N_POS));
R = zeros(size(dt));
Uc = zeros(size(dt));

% Loop through all files
for i = 1:N_POS
    % Read file
    film = readBin('f',i);
    wire = readBin('w',i);
    % Find mean
    Uc(i,:) = mean(wire);
    
    film = lpFilter(film,T,N,50);

    % Adjust for mean and standard deviation
    film = film - mean(film); film = film./std(film);
    wire = wire - mean(wire); wire = wire./std(wire);

    % Pad with zeros
    film_p = [zeros(size(film)); film];
    wire_p = [zeros(size(wire)); wire];
    
    % Compute cross correlation using fft (R_fft goes from -N to N-1)
    R(i,:) = fftshift(ifft(conj(fft(film_p)).*(fft(wire_p))))./N;
end

% Determine dt values to plot
lim = 3*delta./Uc;
dx = -dt.*Uc;
inPlot = abs(dt)<lim;

% Create new grid so that contours can be drawn
[xx, yy] = meshgrid(-3:0.001:3, logspace(log10(y_i), log10(y_f), N_POS));
zz = griddata(dx(inPlot)/delta,y(inPlot),R(inPlot),xx,yy);

figure
contourf(xx,yy/delta,zz);
xlabel('\Deltax/\delta');
ylabel('y/\delta');
axis equal
colorbar

%% Question 5

len = 1000; % Number of values to average either side

% Read and filter values
film = readBin('f',20);
wire = readBin('w',20);
film = lpFilter(film,T,N,50);

% Use two points for gradient
next_point = [film(2:end); 0];
pos_grad = (next_point - film) > 0;

% Average values if gradient of film is positive
condAverage = zeros(2*len+1,1);
nPoints = 0;
for i = len+1:N-len-1
    if (pos_grad(i))
        condAverage = condAverage + wire(i-len:i+len);
        nPoints = nPoints + 1;
    end
end
condAverage = condAverage/nPoints;

plot((-len:len)/Fs,condAverage);
xlabel('\Deltat (s)');
ylabel('Velocity (m/s)');

%% Question 6

% New constants
Fs    = 30*10^3; % Hz
T     = 30;      % s
N     = Fs*T;    % Number of samples
df    = Fs/N;

% Read file
wire = readBin('b',10);
wire = wire - mean(wire);

% Get Fourier transform
g = fft(wire)/N;
g = g(1:N/2+1);
Phi = 2.*g.*conj(g)/df;        % Phi = A^2/2 = 2*g*conj(g)

f = df*(0:N/2)';

% Plot spectrum
figure
semilogx(f,f.*Phi);
xlabel('f (Hz)');
ylabel('f\Phi');

% Check if area under graph equals variance
Area_spectrum = trapz(f,Phi)
variance = mean(wire.^2)

%% Question 7

% New constants
T_total = 30;    % s
Fs    = 30*10^3; % Hz
T     = 7;       % s
N     = Fs*T;    % Number of samples
df    = Fs/N;
FILES = 10;

% Start index of each interval
n_intervals = 1000;
idx = 1:round(Fs*(T_total-T)/n_intervals):Fs*(T_total-T);

% Declare arrays
f   = df*(1:N/2)';
Phi = zeros(N/2,1);
nloops = FILES*n_intervals;
% Loop through all files
for n = 1:FILES
    % Read burst file
    u = readBin('b',n);

    % Loop through intervals with overlapping
    for i = 1:length(idx)
        
        % Extract section of signal
        ui = u(idx(i):N+idx(i)-1);
        ui = ui - mean(ui);
        
        % Get Fourier transform and add contribution
        g   = fft(ui)/N;
        Phi = Phi + 2.*g(2:(N/2+1)).*conj(g(2:(N/2+1)))./(df);
    end
end

% Average power spectral density
Phi = Phi/nloops;

% Plot spectrum
figure
semilogx(f,f.*Phi);
xlabel('f (Hz)');
ylabel('f\Phi');
