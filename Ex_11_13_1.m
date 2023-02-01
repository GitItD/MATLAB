clear
clc

fr = 1000;      % sampling frequency of 1 kHz
sr = 1/fr;      % sampling period
L = 10000;       % if the signal duration is 10 s, then by definition of hz, 1 khz = 10 * 1000 = 10000
t = (0:L-1)*sr;      % time vector = 0:9999/1000 = 0:9.9 == 0:1/1000:10 = 999/10

% The amplitude of our complex sinewave equation
a1 = 0.1;         
a2 = 0.25;
a3 = 0.5;
a4 = 0.75;

% The frequency in hz of our complex sinewave equation
f1 = 1;          
f2 = 10;
f3 = 20;
f4 = 30;

% Our complex sinewave series
csw1 = a1*sin(2*pi*f1*t); 
csw2 = a2*sin(2*pi*f2*t); 
csw3 = a3*sin(2*pi*f3*t);
csw4 = a4*sin(2*pi*f4*t);  

sX1 = csw1 + csw2 + csw3 + csw4 + 2*randn(size(t));     % Add a little bit of noise signal to complex sinewave equation
sX2 = csw1 + csw2 + csw3 + csw4 + 20*randn(size(t));     % Add a lot of noise signal to complex sinewave equation

figure(1), clf
plot(t, sX1)
ylim([0, 35])
title('Small noise time domain')
ylabel('Hz')
xlabel('time (sec)')

figure(2), clf
plot(t, sX2)
ylim([0, 35])
title('Large noise time domain')
ylabel('Hz')
xlabel('time (sec)')

% compute the Fourier transform of the small noise signal
Y1 = fft(sX1);       
P2 = abs(Y1/L);  % the absolute value of Fourier transform divided by signal duration is spectrum 2
P1= P2(1:L/2+1);        % computing spectrum 1 based on P2 and even-valued signal length L
P1(2:end - 1) = 2*P1(2:end - 1); 
% compute the Fourier transform of the large noise signal
Y2 = fft(sX2);       
P4 = abs(Y2/L);  % the absolute value of Fourier transform divided by signal duration is spectrum 4
P3= P4(1:L/2+1);        % computing spectrum 3 based on P4 and even-valued signal length L
P3(2:end - 1) = 2*P3(2:end - 1); 

%{ 
Define the frequency domain f and plot the single-isded amplitude spectrum
P1. The amplitudes are not exactly at 0.7 and 1 as expected, becaue of the
added noise. On average, longer signals produce better frequency
approximations.
%}
figure(3), clf
fs = fr*(0:(L/2))/L;
plot(fs, P1)
title('Single-sided amplitude spectrum of small noise')
xlabel('f (Hz)')
ylabel('|P1(f)|')

figure(4), clf
plot(fs, P3)
title('Single-sided amplitude spectrum of large noise')
xlabel('f (Hz)')
ylabel('|P1(f)|')

%% Answer:

srate  = 1000;
frex   = [ 3  10 15 30 ];
amplit = [ 5  15  5  7 ];
phases = [  pi/8  pi  pi/2  -pi/4 ];
time   = 0:1/srate:10;

% create sine waves
sine_waves = zeros(length(frex),length(time));

% loop over frequencies and define each sine wave
for fi=1:length(frex)
    sine_waves(fi,:) = amplit(fi) * sin(2*pi*time*frex(fi) + phases(fi));
end

% Add (1) a little bit of noise and (2) a lot of noise to make two time series.
% I defined the noise level relative to the average sine wave amplitude
sinesWithSmallNoise = sum(sine_waves,1) + randn(1,length(time)) * mean(amplit);
sinesWithLargeNoise = sum(sine_waves,1) + randn(1,length(time)) * 20*mean(amplit);

% Plot the individual sine waves and the summed time series in separate figures.

figure(5), clf
for fi=1:length(frex)
    subplot(4,1,fi)
    plot(time,sine_waves(fi,:))
    xlabel('Time (s)'), ylabel('Amplitude')
    title([ 'Sine wave of ' num2str(frex(fi)) ' Hz' ])
end

figure(6), clf
subplot(311)
plot(time,sum(sine_waves,1))
xlabel('Time (s)'), ylabel('Amplitude')
title('Average sine wave, NO NOISE')

subplot(312)
plot(time,sinesWithSmallNoise)
xlabel('Time (s)'), ylabel('Amplitude')
title('Average sine wave, LITTLE NOISE')

subplot(313)
plot(time,sinesWithLargeNoise)
xlabel('Time (s)'), ylabel('Amplitude')
title('Average sine wave, LOTS OF NOISE')


labelz = {'no';'little';'lots of'};
figure(7), clf

for i=1:3

    % get no noise, little noise, or lots of noise
    switch i
        case 1
            thing2plot = sum(sine_waves,1);
        case 2
            thing2plot = sinesWithSmallNoise;
        case 3
            thing2plot = sinesWithLargeNoise;
    end

    % plot time domain signal
    subplot(2,3,i)
    plot(time,thing2plot);
    xlabel('Time (ms)')
    title([ 'Time series with ' labelz{i} ' noise' ])

    subplot(2,3,i+3)
    f = fft(thing2plot)/length(time);
    frequencies = linspace(0,srate/2,floor(length(f)/2));
    bar(frequencies,abs(f(1:floor(length(time)/2))).*2)
    set(gca,'xlim',[0 35],'ylim',[0 max(amplit)])
    title([ 'Power spectrum with ' labelz{i} ' noise' ])
    xlabel('Frequencies (Hz)'), ylabel('amp.')
end

%% Answer comparison

fr = 1000;     
sr = 1/fr;      
L = 10000;   
phases = [pi/8, pi, pi/2, -pi/4];       % this was not added in the original
t = (0:L-1)*sr;      % time vector, THIS IS EQUAL TO THE ANSWER

% The amplitude of our complex sinewave equation
a1 = 0.1;         
a2 = 0.25;
a3 = 0.5;
a4 = 0.75;

% The frequency in hz of our complex sinewave equation
f1 = 1;          
f2 = 10;
f3 = 20;
f4 = 30;

% Our complex sinewave series
csw1 = a1*sin(2*pi*f1*t); 
csw2 = a2*sin(2*pi*f2*t); 
csw3 = a3*sin(2*pi*f3*t);
csw4 = a4*sin(2*pi*f4*t);  

sX1 = csw1 + csw2 + csw3 + csw4 + 2*randn(size(t));     % Add a little bit of noise signal to complex sinewave equation
sX2 = csw1 + csw2 + csw3 + csw4 + 20*randn(size(t));     % Add a lot of noise signal to complex sinewave equation

figure(1), clf
plot(t, sX1)
ylim([0, 35])
title('Small noise time domain')
ylabel('Hz')
xlabel('time (sec)')

figure(2), clf
plot(t, sX2)
ylim([0, 35])
title('Large noise time domain')
ylabel('Hz')
xlabel('time (sec)')

% compute the Fourier transform of the small noise signal
Y1 = fft(sX1);       
P2 = abs(Y1/L);  % the absolute value of Fourier transform divided by signal duration is spectrum 2
P1= P2(1:L/2+1);        % computing spectrum 1 based on P2 and even-valued signal length L
P1(2:end - 1) = 2*P1(2:end - 1); 
% compute the Fourier transform of the large noise signal
Y2 = fft(sX2);       
P4 = abs(Y2/L);  % the absolute value of Fourier transform divided by signal duration is spectrum 4
P3= P4(1:L/2+1);        % computing spectrum 3 based on P4 and even-valued signal length L
P3(2:end - 1) = 2*P3(2:end - 1); 

%{ 
Define the frequency domain f and plot the single-isded amplitude spectrum
P1. The amplitudes are not exactly at 0.7 and 1 as expected, becaue of the
added noise. On average, longer signals produce better frequency
approximations.
%}
figure(3), clf
fs = fr*(0:(L/2))/L;
plot(fs, P1)
title('Single-sided amplitude spectrum of small noise')
xlabel('f (Hz)')
ylabel('|P1(f)|')

figure(4), clf
plot(fs, P3)
title('Single-sided amplitude spectrum of large noise')
xlabel('f (Hz)')
ylabel('|P1(f)|')