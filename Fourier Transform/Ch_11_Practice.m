srate = 1000;
frex = [3 10 15 30];
amlit = [ 5 15 5 7];
phases = [pi/8 pi pi/2 -pi/4];
time = 0:1/srate:10;

sine_waves = zeros(length(frex), length(time));

for fi = 1:length(frex)
    sine_waves(fi,:) = amlit(fi) * sin(2*pi*time*frex(fi) + phases(fi));
end

sinesWithSmallNoise = sum(sine_waves, 1) + randn(1, length(time)) * mean(amlit);
sineswithLargeNoise = sum(sine_waves, 1) + randn(1, length(time)) * 20*mean(amlit);

figure(1), clf
for fi = 1:length(frex)
    subplot(4, 1, fi)
    plot(time, sine_waves(fi, :))
    xlabel('Time (s)'), ylabel('Amplitude')
    title(['Sine wave of ' num2str(frex(fi)) ' Hz'])
end

figure(2), clf
subplot(311)
plot(time, sum(sine_waves, 1))
xlabel('Time (s)'), ylabel('Amplitude')
title('Average sine wave, NO NOISE')

subplot(312)
plot(time, sinesWithSmallNoise)
xlabel('Time (s)'), ylabel('Amplitude')
title('Average sine wave, LITTLE NOISE')

subplot(313)
plot(time, sinesWithLargeNoise)
xlabel('Time (s)'), ylabel('Amplitude')
title('Average sine wave, LOTS OF NOISE')

labelz = {'no'; 'little'; 'lots of'};
figure(3), clf

for i = 1:3
    
    switch i
        case 1
            thing2plot = sum(sine_waves, 1);
        case 2
            thing2plot = sinesWithSmallNoise;
        case 3
            thing2plot = sinesWithLargeNoise;
    end
   
    subplot(2, 3, i)
    plot(time, thing2plot);
    xlabel('Time (ms)')
    title(['Time series with ' labelz{i} ' noise'])
    
    subplot(2, 3, i+3)
    f = fft(thing2plot)/length(time);
    frequencies = linspace(0, srate/2, floor(length(f)/2));
    bar(frequencies, abs(f(1:floor(length(time)/2))).*2)
    set(gca, 'xlim', [0 35], 'ylim', [0 max(amplit)])
    title(['Power spectrum with ' labelz{i} ' noise'])
    xlabel('Frequencies (Hz)'), ylabel('amp.')
end