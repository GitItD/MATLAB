%% 23.6.6
%{
Reproduce the left panel of figure 23.6 for other frequencies, ranging
from 2 Hz to 50 Hz. At each iteration of this loop over frequencies,
compute the FFT of this spike-triggered average, and extract power
from the Fourier coefficient closest to the frequency of the wavelet
used to filter the LFP data. Finally, plot the power spectrum as a function
of wavelet frequency. How do these results compare to figure 23.7,
and are you surprised?

2 Hz = 1/2 = 0.5 sec
50 Hz = 1/50 = 0.02 sec

Compute the FFT of spike-average
Extract power from Fourier coefficient closest to frequency o wavelet
%}
cd(['C:\Users\dun81\OneDrive\Documents\Dun II\ET Online\Library Online\Psychology\Computational Neuroscience\' ...
    'MATLAB for Brain and Cognitive Scientists - The MIT Press (2017)\MATLAB_BrainAndCognitiveScientists-main\ch23\'])
load spikefieldData.mat

% Parameters
wherespikes = find(spikeTimes);
win = 400;
wherespikes(wherespikes<win | wherespikes>length(lfp)-win) = [];
spikeLFP = zeros(length(wherespikes), win*2+1);
timevec = (-win:win)/(1000/srate);

% Find spikes
for si = 1:length(wherespikes)
    spikeLFP(si,:) = lfp(wherespikes(si) - win : wherespikes(si) + win);
end

%{
At each iteration of this loop over frequencies,
compute the FFT of this spike-triggered average, and extract power
from the Fourier coefficient closest to the frequency of the wavelet
used to filter the LFP data.
%}

for ii = 1:length(timevec)
    avgspikeLFP = mean(spikeLFP);
    FFTspikeLFP = fft(avgspikeLFP);
    FFTpower(ii,:) = abs((FFTspikeLFP(floor(1:length(timevec)/2) + 1)).^2);
end

%{ 
Create wavelet
csine = exp(2*1i*pi*freq*spikeLFP);
gaus = exp( -(spikeLFP.^2) / (2*s^2) );
cmw = csine .* gaus;
[~, idx] = min(abs(freqs - waveletFreq));
%}

%FFTWav = (bsxfun(@times, FFTspikeLFP, sin(pi*win*spikeLFP)).^2);


%FFTWav = bsxfun(@rdivide, FFTWav, avgspikeLFP);

figure(1), clf
% subplot of FFTAvg
subplot(121)
timevec = (win-380:win)/(1000/srate);
h = plot(timevec, spikeLFP(1:length(timevec),:));
hold on
plot(timevec, (avgspikeLFP(1:length(timevec))), 'k', 'linew', 4)
set(gca, 'xlim', timevec([1, end]))
set(h, 'color', [.7, .7, .7])
xlabel('Peri-spike time(ms)'), ylabel('Trials')
% subplot of Wavelet Power Spectrum
subplot(122)
plot(timevec, FFTpower(1:length(timevec)))
xlabel('Frequencies (Hz)'), ylabel('Amplitude')
title('Power spectrum of LFP filter')
xlim([timevec(end) - length(timevec), length(timevec)]);

%% 23.6.6 ChatGPT answer
freqs = 2:50; % range of frequencies
waveletWidth = 7; % width of each wavelet in cycles
spikeAvg = mean(spikeLFP);
% initialize output matrix
waveletPower = zeros(length(freqs),1);

% Compute the FFT of the spike-triggered average
spikeAvgFFT = fft(spikeAvg);

% Find the frequency resolution of the FFT
freqRes = srate / length(spikeAvg);

% Compute the frequencies corresponding to each FFT coefficient
freqsFFT = 0:freqRes:srate-freqRes;

% Loop over frequencies
for fi = 1:length(freqs)
    % Create a wavelet for the current frequency
    waveletFreq = freqs(fi);
    nCycles = waveletWidth*waveletFreq;
    time = -nCycles/2:nCycles/2;
    wavelet = exp(2*1i*pi*waveletFreq*time/srate) .* exp(-time.^2/(2*(waveletWidth*waveletFreq/(2*pi*srate))^2));
    
    % Convolve the wavelet with the spike-triggered average
    spikeLFP = conv(spikeAvg,wavelet,'same');
    
    % Compute the FFT of the convolved signal
    spikeLFPFFT = fft(spikeLFP);
    
    % Find the index of the FFT coefficient closest to the frequency of the wavelet
    [~, idx] = min(abs(freqsFFT - waveletFreq));
    
    % Compute the power of the FFT coefficient closest to the frequency of the wavelet
    waveletPower(fi) = abs(spikeLFPFFT(idx))^2;
end

% Plot the power spectrum as a function of wavelet frequency
figure(2), clf
plot(freqs, waveletPower,'o-')
xlabel('Frequency (Hz)'), ylabel('Power')

%% 23.6.7
%{
Here is another way to compute spike-field coherence: discretize the
LFP phase-angle time series into N bins (e.g., 10 bins), and compute
the average number of spikes per bin. Then plot the results as a bar
plot of spike count (y axis) as a function of LFP phase (x axis). Try this
for a few frequencies and confirm that the results are consistent with
figure 23.7.
%}


%% 23.6.8
%{
The spike-field coherence value at 8 Hz is around 0.32 (the exact value
will depend on parameters). Can this be considered statistically significant?
Generate a null hypothesis distribution by randomizing spike
times and recomputing phase clustering. Perform 1,000 iterations, plot
the null hypothesis distribution and the observed value, and compute
a p-value. This method is actually too liberal, because the temporal
characteristics of the shuffled action potentials are no longer physiologically
plausible. Therefore, also try the Ågcut-and-shiftÅh procedure, in
which the phase-angle time series is cut in a random location on each
iteration, and the second half is placed before the first half. How do the
two null-hypothesis distributions compare?
%}