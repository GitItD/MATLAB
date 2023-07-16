%% Matlab for Brain and Cognitive Scientists
% Matlab code for Chapter 18
% Mike X Cohen
% 
% This code accompanies the book 
%      "Matlab for Brain and Cognitive Scientists" (MIT Press, 2017).
% Using the code without following the book may lead to confusion, 
% incorrect data analyses, and misinterpretations of results. 
% Mike X Cohen assumes no responsibility for incorrect use of this code. 

%% frequency resolution
%{
The frequency resolution is a measure of the smallest frequency difference
that can be distinguished by the analysis. It is an important consideration in 
the analysis of signals in the frequency domain, as it determines the accuracy 
and precision of the frequency content of the signal.
%}

srate = 739;  % sampling rate in Hz
n     = 1234; % number of time points

%{
This line of code generates a vector hz containing a range of frequencies from 0 Hz to half of the sampling rate, 
with a frequency resolution that depends on the number of time points in
the following format:

hz = linspace(start_time, end_time, Nyquist_frequency)
%}
hz = linspace(0,srate/2,floor(n/2)+1);

%{
Frequency resolution can be measured in different ways depending on the context and application. 
In some cases, a specific method may be more appropriate or useful than others. 
For example, the average derivative of frequencies may be a more appropriate measure 
when the frequency range is large and the frequency resolution is variable, 
while the difference between any two frequencies may be more useful when a specific 
frequency range is being studied. The ratio of sampling rate to number of data points is a 
simple measure that can be used when a rough estimate of frequency resolution is sufficient. 
Ultimately, the choice of method depends on the specific context and goals of the analysis.
%}
% frequency resolution, measured in several different ways:
% average derivative of frequencies in Hz
freqres = mean(diff(hz));

% difference between any two frequencies
freqres = hz(2)-hz(1);

% ratio of sampling rate to number of data points
freqres = srate/n;

%% finding the right N given a desired frequency resolution
%{
This code snippet calculates the number of time points required (nFFT) to 
achieve a desired frequency resolution (freqres) given a sampling rate (srate).

The formula used is nFFT = ceil(srate/freqres), which calculates the number of 
time points needed to achieve a frequency resolution of freqres using the given sampling rate. 
The ceil() function is used to round up to the nearest integer, as the number of time points must be an integer.

For example, if the sampling rate is 1000 Hz and the desired frequency resolution is 0.25 Hz, 
the number of time points needed to achieve this frequency 
resolution is 4000 (nFFT = ceil(1000/0.25) = 4000).
%}
srate   = 1000;
freqres = .25; % in Hz

% notice that this line is the same as (freqres=srate/n) but solved for n
nFFT    = ceil( srate/freqres );

%% difference between sampling rate and number of time points for Fourier frequencies
%{
This code cell shows the difference between frequency resolution and the Nyquist frequency. 
It generates three signals with different sampling rates and time durations, and computes their Fourier transforms.

In the time domain, each signal is displayed as a waveform. In the frequency domain, the amplitude of each frequency 
component is plotted against its frequency, up to the Nyquist frequency, which is half of the sampling rate. 
The Nyquist frequency is important because it is the highest frequency that can be represented in the signal, 
and any frequency above this will be aliased, meaning that it will appear as a lower frequency.

This code cell illustrates the difference between frequency resolution and the Nyquist frequency by showing
 how changing the sampling rate and the number of time points affect the resolution and the maximum 
frequency that can be represented in the signal.
%}
% This cell illustrates the difference between 
% frequency resolution and the Nyquist frequency.

% temporal parameters
srates  = [100 100 1000];
timedur = [  1  10    1];

% a few more specifications
freq    =     5; % in Hz
colors  = 'kmb';
symbols = 'op.';


figure(1), clf
legendText = cell(size(timedur));
for parami=1:length(colors)
    
    % define sampling rate in this round
    srate = srates(parami); % in Hz
    
    % define time
    time = -1:1/srate:timedur(parami);
    
    % create signal (Morlet wavelet)
    signal = cos(2*pi*freq.*time) .* exp( (-time.^2) / .05 );
    
    % compute FFT and normalize
    signalX = fft(signal)/length(signal);
    signalX = signalX./max(signalX);
    
    % define vector of frequencies in Hz
    hz = linspace(0,srate/2,floor(length(signal)/2)+1);
    
    
    % plot time-domain signal
    subplot(211)
    plot(time,signal,[colors(parami) symbols(parami) '-'],'markersize',10,'markerface',colors(parami)), hold on
    set(gca,'xlim',[-1 1])
    xlabel('Time (s)'), ylabel('Amplitude')
    title('Time domain')
    
    % plot frequency-domain signal
    subplot(212), hold on
    plot(hz,abs(signalX(1:length(hz))),[colors(parami) symbols(parami) '-'],'markersize',10,'markerface',colors(parami))
    xlabel('Frequency (Hz)'), ylabel('Amplitude')
    title('Frequency domain')
    
    legendText{parami} = [ 'srate=' num2str(srates(parami)) ', N=' num2str(timedur(parami)+1) 's' ];
end

legend(legendText)
zoom on

%% edges in the time domain produce "artifacts" in the frequency domain
%{
In the above example, a time series signal is created with a sharp edge by setting a range 
of time points to a high value of 10. The plot in the upper panel shows the time-domain 
signal with the sharp edge. The plot in the lower panel shows the power spectrum of the 
signal obtained by taking the absolute value of the FFT (Fast Fourier Transform) of the signal, 
squaring the amplitude spectrum and normalizing it by dividing by the number of data points 
(in this case, 100). The x-axis of the lower plot shows the frequency in arbitrary units, and 
the y-axis shows the power in arbitrary units.

In this particular case, the power spectrum shows a clear peak at a frequency of 0.5, 
which corresponds to the inverse of the distance between the two edges of the signal 
(100 data points - 5 data points with high amplitude = 95 data points, and the inverse 
of 95 is approximately 0.5). This illustrates the relationship between the sharpness of 
the edges in the time domain and the spread of power in the frequency domain.
%}
% create time series with sharp edge
ts = zeros(100,1);
ts(48:52) = 10;


figure(2), clf

% plot time-domain signal
subplot(211)
plot(ts,'k','linew',3)
xlabel('Time (a.u.)'), ylabel('Amplitude (a.u.)')
set(gca,'ylim',[-.5 10.5])

% plot power spectrum of signal
subplot(212)
plot(abs(fft(ts)/100).^2,'ro','linew',1.5,'markerface','g')
xlabel('Frequencies (a.u.)'), ylabel('Power (a.u.)')
%{
The fft function is used to compute the Fourier transform of the time series, 
which yields a complex-valued output. The abs function is then used to compute 
the magnitude of the Fourier coefficients, and the .^2 operator is used to compute 
their squared values, which represent the power at each frequency.

Since the time series has a sharp edge, its Fourier transform exhibits high-frequency 
components that correspond to the abrupt transition from zero to the non-zero value. 
These high-frequency components are represented by the peaks in the power spectrum 
that appear around frequencies 0.5 and -0.5 (since the time series is real-valued, 
its Fourier transform is symmetric and its frequencies are centered around zero
%}

%% illustration of tapering
%{
Tapering is a way to attenuate the edges of the data to reduce spectral leakage. 
Spectral leakage occurs because the Fourier transform assumes that the signal is 
periodic over the duration of the time series, which is not always the case. 
When the signal is not periodic, the Fourier transform will cause the energy to 
be spread across multiple frequency bins, making it difficult to identify the true 
spectral content of the signal. Tapering can be achieved by multiplying the data with 
a window function that gradually attenuates the edges of the data. In this figure, the 
original random data is plotted in blue, the tapered data is plotted in red, and the Hann window is 
plotted in black. As can be seen, the tapered data has attenuated edges compared to the 
original data, which is consistent with the Hann window.
%}
figure(3), clf

n = 100;
r = randn(1,n);
hannwin = .5*( 1-cos(2*pi*(0:n-1)/(n-1)) );
plot(r), hold on
plot(r.*hannwin,'r')
plot(hannwin,'k--')
%{
The variable hannwin is a vector that defines a Hann (Hanning) window of length n. 
A Hann window is a type of tapering function that is often used in signal processing 
to smooth the edges of a signal or to reduce spectral leakage in frequency analysis. 
The Hann window is defined as:

hannwin(n) = 0.5 * (1 - cos(2 * pi * (0:n-1) / (n-1)))

where n is the length of the window. The Hann window tapers smoothly from 0 at 
both ends to 1 in the middle, and has a bell-shaped profile that is similar to a Gaussian 
distribution. By multiplying a signal with a Hann window of the same length, the signal 
is smoothly tapered to zero at both ends, which can help reduce the effects of spectral 
leakage and improve frequency resolution in Fourier analysis.
%}

xlabel('Time (a.u.)')
ylabel('Amplitude (a.u.)')
title('Tapering tapers the data')

legend({'original data';'tapered data';'Hann taper'})

%% epoching continuous data
%{
Epoching is a common data processing step in many fields, particularly in neuroscience, 
where it is used to extract meaningful time series data from continuous recordings. 
In essence, epoching involves dividing a long continuous signal into shorter segments or 
"epochs," typically with a fixed duration, to facilitate analysis and interpretation of the data.

Epoching is particularly useful in situations where the signal of interest occurs intermittently 
or is otherwise hard to distinguish from noise or other artifacts. By breaking the signal into 
smaller segments, researchers can more easily identify and extract the relevant information, 
and discard any portions of the signal that are not relevant to the analysis.

This code generates a continuous time series of randomly generated data (data) with a 
sampling rate of 512 Hz and a length of 21 seconds. It also specifies a desired epoch length 
(epochLms) of 2345 ms and converts it to the corresponding number of data points (epochLidx).

It then calculates the number of epochs that can be created from the available data, 
given the desired epoch length, and generates a matrix (epochs) in which each row 
represents an epoch of length epochLidx. If the number of data points is not evenly 
divisible by epochLidx, the last epoch will contain only the remaining data points.

This code is useful for creating epochs of data for use in analyses such as 
event-related potentials (ERPs) or time-frequency analyses.
%}
srate = 512; % in hz
n     = 21*srate; % 21 seconds of â€œactivityâ€?

% continuous data
data  = randn(1,n);

% desired length of epoch
epochLms  = 2345; % in ms
% and convert to indices (time units)
epochLidx = round( epochLms / (1000/srate) );

% how many epochs will we have?
nEpochs = floor(n/epochLidx);

% finally, create the epochs using however much data is possible
epochs  = reshape( data(1:nEpochs*epochLidx) ,nEpochs,epochLidx);

%% using the reshape function
%{
The reshape function in MATLAB is used to change the shape (i.e., the number of rows and columns) 
of an array, without changing the total number of elements in the array.

The general syntax for reshape is:
B = reshape(A, [m, n]);
where A is the input array, m and n are the desired number of rows and columns in the output array, 
respectively, and B is the reshaped output array.
%}

% original data
data = 1:12;

% various reshaping options
reshape(data,1,[])
reshape(data,[],1)
reshape(data,3,4)
reshape(data,4,3)
reshape(data,4,4)

%% reshaping problems

nchans = 10;
data = randn(nchans,n);
epochs = reshape(data(:,1:nEpochs*epochLidx),[nchans nEpochs epochLidx]);


figure(4), clf

subplot(211)
plot(data(1,1:100),'r'), hold on
plot(squeeze(epochs(1,1,1:100)))
xlabel('Time (a.u.)'), ylabel('Amplitude (a.u.)')
title('Uh oh.')

subplot(212)
epochs1 = reshape(data(:,1:nEpochs*epochLidx),[nchans epochLidx nEpochs]);
plot(data(1,1:100),'r','linew',2), hold on
plot(squeeze(epochs1(1,1:100,1)),'ko--')
xlabel('Time (a.u.)'), ylabel('Amplitude (a.u.)')
title('Looks better.')

%% finding the right frequency
%{

%}
srate = 1000;
nFFT  = 5128;
hz    = linspace(0,srate/2,floor(nFFT/2)+1);

% does this give frequencies between 8 Hz and 12 Hz?
hz([8 12])


figure(5), clf

desfreq = 8; % desired frequency in Hz

subplot(131)
plot(hz), hold on
plot(get(gca,'xlim'),[desfreq desfreq],'k--')
xlabel('Indices into variable hz')
ylabel('Actual frequencies in Hertz')
set(gca,'xlim',[0 400],'ylim',[0 30])
title('hz')
axis square


subplot(132)
plot(hz-desfreq), hold on
plot(get(gca,'xlim'),[0 0],'k--')
xlabel('Indices into variable hz')
ylabel([ 'Hertz relative to ' num2str(desfreq) ])
set(gca,'xlim',[0 400],'ylim',[-30 30])
title('hz-desfreq')
axis square

subplot(133)
plot(abs(hz-desfreq)), hold on
plot(get(gca,'xlim'),[0 0],'k--')
xlabel('Indices into variable hz')
ylabel([ 'Hertz relative to ' num2str(desfreq) ])
set(gca,'xlim',[0 400],'ylim',[-30 30])
title('abs(hz-desfreq)')
axis square

%% power spectrum of real EEG data before and after tapering
%{
The code loads EEG data (one electrode from a human during resting-state) and computes its power spectrum. 
It then plots the power spectrum before and after applying a Hann taper to the data, which tapers the edges of 
the data to reduce spectral leakage in the Fourier transform.

The first subplot in the figure shows the power spectrum of the raw time series in black and the power 
spectrum of the tapered time series in red. The second subplot shows the normalized power spectra 
of the raw and tapered time series, with the same color coding. The x-axis represents frequency in Hz, 
and the y-axis represents power in microvolts squared or normalized power.
%}
% load EEG data (one electrode from a human during resting-state)
load EEGrestingState
npnts = size(eegdata,1);

dataX = fft(eegdata,[],1)/npnts;
hz = linspace(0,srate/2,floor(npnts/2)+1);

figure(5), clf
subplot(211)
pow = mean( abs(dataX(1:length(hz),:)).^2 ,2);
plot(hz,pow,'k.-')
%{
In the code pow = mean( abs(dataX(1:length(hz),:)).^2 ,2);, dataX is the 
Fourier transform of the EEG data eegdata after applying the tapering window. 
The length(hz) ensures that only the non-negative frequencies are included in the 
power spectrum calculation, and abs(dataX(1:length(hz),:)).^2 calculates the squared 
magnitude of the Fourier coefficients at these frequencies for each time point. 
Finally, mean(...,2) calculates the mean across time points for each frequency, 
resulting in a power spectrum estimate.
%}
subplot(212)
plot(hz,pow./max(pow),'k.-')



% uh oh, no tapering. Better redo:
hannwin = .5 - cos(2*pi*linspace(0,1,npnts))/2;
dataTaper = bsxfun(@times,eegdata,hannwin');
dataX = fft(dataTaper,[],1)/npnts;

subplot(211), hold on
pow = mean( abs(dataX(1:length(hz),:)).^2 ,2);
plot(hz,pow,'r.-')

xlabel('Frequency (Hz)'), ylabel('Power (\muV^2)')
set(gca,'xlim',[0 55])
legend({'Raw time series';'Tapered time series'})

subplot(212), hold on
plot(hz,pow./max(pow),'r.-')

xlabel('Frequency (Hz)'), ylabel('Power (norm.)')
set(gca,'xlim',[0 55])
legend({'Raw time series';'Tapered time series'})

%% extracting power from a specified frequency band
%{
To compute the power within frequency boundaries, you can first use the dsearchn 
function to find the indices in the frequency vector hz that correspond to the desired 
frequency range. Then, you can extract the power values from the Fourier spectrum 
within these indices and average them to get the average power within the frequency range.
%}
desfreq = [8 12]; % frequency boundaries for averaging
idx = dsearchn(hz',desfreq'); % find indices in hz that correspond to desfreq
pow = mean( abs(dataX(idx(1):idx(2))).^2 ); % extract and average power within frequency range

%% temporal non-stationarities and the FFT
%{
What are nonstationary in the time-domain signal?

Nonstationary signals are signals whose properties vary over time. In other words, nonstationary signals are 
time-varying signals. These signals can have characteristics such as amplitude modulation, frequency modulation, 
and time-varying phase. Nonstationary signals can arise from a variety of sources, such as the changing properties 
of a physical system over time or changes in the environment in which the signal is recorded. Examples of nonstationary 
signals include speech signals, electroencephalography (EEG) signals, and music signals.

How does performing time-frequency analyses help with this problem?

Performing time-frequency analyses can help with the problem of nonstationarity in the time-domain signal because 
it allows for the examination of changes in the frequency content of the signal over time. By decomposing the signal into 
its time-varying frequency components, it is possible to identify when and where changes in frequency occur, which can 
provide insight into changes in the underlying neural activity or physiological processes that are generating the signal. 
This can be particularly useful in fields such as neuroscience and biomedical signal processing, where nonstationarity is 
common and changes in frequency content are often of interest.

In the code below, two sine waves are created with the same frequency of 3 Hz but with different amplitudes. 
The amplitude of one sine wave (signal1) linearly increases over time, while the amplitude of the other sine wave 
(signal2) is stationary. The signals are then Fourier transformed, and the power spectrum of each signal is plotted. 
The plot shows that the power spectrum of signal1 varies across frequencies, while the power spectrum of signal2 
is constant across frequencies. The code demonstrates the concept of non-stationarity in signals and how it affects 
their power spectra.
%}
% define some signal parameters
srate = 1000;
t = 0:1/srate:10;
n = length(t);
f = 3; % frequency in Hz

% create time-increasing amplitude, and stationary amplitude
ampl1 = linspace(1,10,n);
% the next line creates a randomly varying amplitude, FYI
% ampl1 = abs(interp1(linspace(t(1),t(end),10),10*rand(1,10),t,'spline'));
ampl2 = mean(ampl1);

% and create the sine waves
signal1 = ampl1 .* sin(2*pi*f*t);
signal2 = ampl2 .* sin(2*pi*f*t);

% FFTs etc.
signal1X = fft(signal1)/n;
signal2X = fft(signal2)/n;
hz = linspace(0,srate/2,floor(n/2)+1);


figure(6), clf
subplot(231)
h(1) = plot(t,signal2);
hold on
plot(t,signal1,'k')
xlabel('Time'), ylabel('amplitude')

subplot(234)
h(2) = plot(hz,2*abs(signal2X(1:length(hz))),'o-');
hold on
plot(hz,2*abs(signal1X(1:length(hz))),'ks-','markerface','k')

xlabel('Frequency (Hz)'), ylabel('amplitude')
set(gca,'xlim',[1 7])
legend({'Stationary';'Non-stationary'})

%% frequency non-stationarity
%{
In this section, the code demonstrates frequency non-stationarity. 
Two sine waves with frequencies 2 Hz and 10 Hz are created. 
Then, the frequency of the second signal is linearly increased from 
2 Hz to 10 Hz across time to create a chirp signal. The stationary and chirp 
signals are then plotted in the time and frequency domains. In the time domain, 
the stationary signal has a constant frequency, while the frequency of the chirp signal linearly 
increases over time. In the frequency domain, the amplitude of the stationary signal is constant across 
all frequencies, while the amplitude of the chirp signal varies depending on the frequency.
%}
% create chirp and frequency-stationary signal
f  = [2 10];
ff = linspace(f(1),mean(f),n);
signal1 = sin(2*pi.*mean(ff).*t);
signal2 = sin(2*pi.*ff.*t);

% FFTs, etc.
signal1X = fft(signal1)/n;
signal2X = fft(signal2)/n;
hz = linspace(0,srate/2,floor(n/2));

% and plot
subplot(232)
h(3) = plot(t,signal2);
hold on
plot(t,signal1,'k')
xlabel('Time'), ylabel('amplitude')
set(gca,'ylim',[-1.1 1.1])

subplot(235)
plot(hz,2*abs(signal2X(1:length(hz))),'k.-'), hold on
h(4) = plot(hz,2*abs(signal1X(1:length(hz))),'.-');
xlabel('Frequency (Hz)'), ylabel('amplitude')
set(gca,'xlim',[0 20])


% Gray lines for a cloudy day. Feel free to adjust 
% the color to match your local weather.
set(h,'color',.6*ones(3,1));

%% even weirder stationarities

% triangle-ish frequency modulation
freqTS = abs(mod(t,2)-1)*10;
meanF = mean(freqTS);
k = 2*pi/srate;
y = sin(2*pi.*meanF.*t + k*cumsum(freqTS-meanF));


subplot(233), plot(t,y)
p = abs(fft(y)/n)*2;
set(gca,'ylim',[-1.1 1.1])

subplot(236), plot(hz,p(1:length(hz)))
set(gca,'xlim',[0 20])

%%

srate = 1000; 
time  = 0:1/srate:5;

% create a signal at 10 Hz with amplitude modulated at 1 Hz
sig = sin(2*pi*time) .* sin(2*pi*10*time);

clf

% plot the time-domain signal
subplot(211)
plot(time,sig)
title('Time domain')
xlabel('Time (s)'), ylabel('amplitude')

% and plot it in the frequency domain
subplot(212)
plot(linspace(0,srate,length(time)), abs(fft(sig)/length(sig)).^2)
set(gca,'xlim',[5 15])
title('Frequency domain')
xlabel('Frequency (Hz)'), ylabel('Power')

% questions: Why don't you see the 10 Hz component in the power spectrum?
%            What happens if you add 1 to the amplitude-modulating function (use parentheses)?

%% spectral coherence
%{
This code generates signals and calculates their spectral coherence. The parameters of the signals are defined, 
including the sampling rate, time vector, number of samples, and frequency vector. Three base signals are also defined, 
each containing random Gaussian noise. Two new signals are then created by adding sine waves with different frequencies 
to the base signals. The spectral coherence between the two new signals is then calculated using the Fourier transform of 
each signal, multiplying the Fourier coefficients of the two signals and taking the absolute value squared. 
The resulting coherence values are then plotted in a graph, along with the time-domain signals and power 
spectra of the original signals.
%}
% define parameters
srate = 1000;
t     = 0:1/srate:9;
n     = length(t);
hz    = linspace(0,srate/2,floor(n/2));

% define frequencies
f = [10 14 8];

% define 'base' signals
f_ts1 = (2*pi*cumsum(5*randn(1,n)))/srate;
f_ts2 = (2*pi*cumsum(5*randn(1,n)))/srate;
f_ts3 = (2*pi*cumsum(5*randn(1,n)))/srate;

% notice how the different signals contain different
% parts of the 'base' signals
sigA = sin(2*pi.*f(1).*t + f_ts1) + randn(size(t));
sigB = sigA + sin(2*pi.*f(2).*t + f_ts2);
sigA = sigA + sin(2*pi.*f(3).*t + f_ts3);
% sigB = sigA;

% spectral coherence
sigAx = fft(sigA)/n;
sigBx = fft(sigB)/n;
specX = abs(sigAx.*conj(sigBx)).^2;
% specX = specX./( abs(sigAx).^2 .* abs(sigBx).^2 );

%%% Note: The normalization is not valid when there is only a single trial. 
%%%  The data should be epoched, or you can use the function mscohere, which
%%%  will epoch continuous data (also try mscohere forcing a single window to
%%%  confirm that the normalized coherence values are all 1's).


figure(7), clf

% plot time-domain signals
subplot(231)
plot(t,sigA)
set(gca,'xlim',t([1 end]))
xlabel('Time (s)'), ylabel('Amplitude')
title('Time-domain signal A')

subplot(232)
plot(t,sigB)
set(gca,'xlim',t([1 end]))
xlabel('Time (s)'), ylabel('Amplitude')
title('Time-domain signal B')


% then plot power spectra
subplot(234)
plot(hz,abs(sigAx(1:length(hz))).*2)
set(gca,'xlim',[0 max(f)*1.3])
xlabel('Frequency (Hz)'), ylabel('Amplitude')
title('Signal A')

subplot(235)
plot(hz,abs(sigBx(1:length(hz))).*2)
set(gca,'xlim',[0 max(f)*1.3])
xlabel('Frequency (Hz)'), ylabel('Amplitude')
title('Signal B')

% finally, plot coherence
subplot(236)
plot(hz,abs(specX(1:length(hz))))
set(gca,'xlim',[0 max(f)*1.3])
xlabel('Frequency (Hz)'), ylabel('Coherence')
title('Spectral coherence: A-B')

%% end.

