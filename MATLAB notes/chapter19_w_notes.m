%% Matlab for Brain and Cognitive Scientists
% Matlab code for Chapter 19
% Mike X Cohen
% 
% This code accompanies the book 
%      "Matlab for Brain and Cognitive Scientists" (MIT Press, 2017).
% Using the code without following the book may lead to confusion, 
% incorrect data analyses, and misinterpretations of results. 
% Mike X Cohen assumes no responsibility for incorrect use of this code. 

%% complex Morlet wavelet
%{
This MATLAB code defines and plots a Complex Morlet wavelet in the time
domain and its power spectrum.

The first section of the code defines the sampling rate, the time vector, the 
frequency of the wavelet, and the parameters for the Gaussian window. 
The Complex Morlet wavelet is then created by multiplying a complex sine 
wave with a Gaussian window in the time domain.

The second section of the code calculates the power spectrum of the 
Morlet wavelet using the FFT. The power spectrum is plotted as a function 
of frequency. The result shows that the Morlet wavelet has a broad 
frequency content centered around the specified frequency 
(frex = 6.5 in this example), which is a characteristic of wavelets 
with a long duration and good time-frequency localization.
%}
srate = 1000;
wavtime = -2:1/srate:2;
frex = 6.5;
s = 5/(2*pi*frex);
csine = exp(2*1i*pi*frex*wavtime);
gaus = exp( -(wavtime.^2) / (2*s^2) );
cmw = csine .* gaus;

figure(1), clf
plot3(wavtime, real(cmw), imag(cmw))
xlabel('Time'), ylabel('Real part'), zlabel('Imaginary part')
axis image
title('Complex Morlet wavelet in the time domain')
rotate3d on


figure(2), clf
hz = linspace(0, srate/2, floor(length(wavtime)/2)+1);
cmwX = fft(cmw)/length(wavtime);
plot(hz, abs(cmwX(1:length(hz)))*2)
set(gca, 'xlim', [0 frex*3])
xlabel('Frequencies (Hz)'), ylabel('Amplitude')
title('Power spectrum of Morlet wavelet')

%{
The code you provided computes the power spectrum of the complex Morlet wavelet defined earlier. 
It first defines a frequency vector hz using linspace that ranges from 0 to the Nyquist frequency 
(half the sampling rate) and has the same number of points as half the length of the time-domain wavelet.

The fft function is then used to compute the Fourier transform of the complex Morlet wavelet in the time domain, 
which is then normalized by dividing by the length of the time-domain signal.

The abs function is used to compute the absolute value of the Fourier coefficients, 
and multiplying by 2 gives the power spectrum. Finally, the plot is created using the plot function, 
with the x-axis set to the frequency vector hz and the y-axis set to the power spectrum. 
The xlim function is used to set the x-axis limits to show only frequencies up to three times 
the center frequency of the wavelet.

The resulting plot shows the frequency-domain representation of the Morlet wavelet, 
which is a bell-shaped curve centered at the wavelet's central frequency (frex in the code).
%}

%% Convolution!
%{
This code is performing a time-frequency analysis on a linear chirp signal using convolution.
In signal processing, a chirp is a signal in which the frequency increases or decreases over time. 
It can be represented as a sinusoidal waveform whose frequency varies with time. 
A linear chirp is a signal whose frequency increases or decreases linearly with time, 
while a nonlinear chirp has a more complex frequency modulation.

First, the signal parameters are defined. sigtime is a time vector from 0 to 10 seconds 
with a sampling rate of srate. A linear chirp signal is created using sine waves with increasing
 frequency from 2 Hz to 10 Hz.

Next, the convolution parameters are defined. nData is the length of the signal, nKern is 
the length of the wavelet time vector wavtime, and nConv is the length of the resulting 
convolution. nHfkn is half of the length of the wavelet time vector plus one.

The convolution is performed by taking the inverse Fourier transform of the product 
of the Fourier transforms of the continuous Morlet wavelet cmw and the signal signal. 
The resulting time series is then truncated by removing the edge effects.

Finally, two subplots are created. The first subplot shows the original signal in the time 
domain, while the second subplot shows the time-frequency power of the signal at a 
single frequency defined by frex.
%}
% some signal parameters
sigtime = 0:1/srate:10;
n = length(sigtime);

% create a linear chirp
f  = [2 10];
ff = linspace(f(1),mean(f),n);
signal = sin(2*pi.*ff.*sigtime);
%{
This code creates a linear chirp signal using the following parameters:

f: A 1x2 vector defining the start and end frequencies of the chirp signal in Hz.
n: The number of time points for the signal.
sigtime: A 1xn vector defining the time points for the signal in seconds.
The code first generates a linearly spaced vector of frequencies (ff) between the 
start frequency and the mean of the start and end frequencies. Then it generates a 
sinusoidal signal at each time point with a frequency determined by ff. 
This results in a signal that starts at a low frequency and linearly increases 
to a higher frequency.
%}

% convolution parameters
nData = length(sigtime);
nKern = length(wavtime);
nConv = nData+nKern-1;
nHfkn = floor(length(wavtime)/2)+1;

% the meat and potatoes of convolution
as = ifft( fft(cmw,nConv) .* fft(signal,nConv) );
as = as(nHfkn:end-nHfkn+1);
%{
This code is performing the convolution between the signal and the cmw using the Fourier transform.

The fft function is used to convert the signal and the convolution kernel cmw into the frequency domain. 
The multiplication of the two spectra in the frequency domain is equivalent to the convolution of the two 
signals in the time domain. The ifft function is then used to convert the convolution result back to the time domain.

The nConv variable is used to determine the length of the convolution result, which is the sum of the length 
of the signal and the length of the kernel minus 1. The nHfkn variable is used to determine the center of the 
kernel to be used in the final output.

Finally, the as variable contains the convolution result after discarding the edges resulting from the convolution 
process using nHfkn.
%}

% plot the time-domain signal (chirp)
figure(3), clf
subplot(211)
plot(sigtime, signal)
set(gca,'ylim',[-1.1 1.1])
xlabel('Time (s)'), ylabel('Amplitude')
title('Time-domain signal')

% and the power time series at one frequency
subplot(212)
plot(sigtime, abs(as).^2)
xlabel('Time (s)'), ylabel('Amplitude')
title([ 'Time-frequency power (' num2str(frex) ' Hz)' ])

%% now a time-frequency plane
%{
This code generates a time-frequency representation of a signal using Morlet wavelets. 
It initializes a time-frequency matrix tf, and for each frequency value in frex, it creates a 
corresponding Morlet wavelet using the parameters s and wavtime. It then computes 
the FFT of the signal (sigX) and the FFT of the wavelet (cmwX), multiplies them, and takes 
the IFFT to obtain the convolution of the signal with the wavelet. The absolute value of this 
convolution is stored in the tf matrix for the corresponding frequency and time points. 
Finally, the contourf function is used to plot the time-frequency representation of the signal.

In this example, the Morlet wavelet was created using a loop and simple mathematical 
operations to define the wavelet in the time-domain. In previous example, the Morlet wavelet is 
defined using a mathematical formula that includes a Gaussian function and a complex 
exponential. Previous formula is used to directly create the wavelet in the time-domain without 
the need for a loop and its approach is more efficient and faster, especially when creating 
many wavelets of different frequencies and scales.
%}
% new frequency parameters
nFrex = 30;
frex = linspace(1,15,nFrex);

% width of Gaussian
s = linspace(4,12,nFrex) ./ (2*pi.*frex);

% initialize output matrix
tf = zeros(nFrex, length(signal));

for fi=1:nFrex
    
    % Take the FFT of the signal.
    % Does it make sense to have this code *inside* the frequency loop?
    sigX = fft(signal, nConv);
    
    % create Morlet wavelet
    cmw = exp(  2*1i*pi*frex(fi)*wavtime + ...
              -(wavtime.^2)/(2*s(fi)^2) );
    
    % compute its FFT      
    cmwX = fft(cmw,nConv);
    
    % 'meat' of convolution
    as = ifft( sigX .* cmwX );
    tf(fi,:) = abs( as(nHfkn:end-nHfkn+1) )*2;
end

% and plot
figure(4), clf
contourf(sigtime,frex,tf,40,'linecolor','none')

%% check the frequency representation of the wavelets
%{
This code creates a plot of the wavelet in the time domain and its power spectrum 
in the frequency domain for each frequency in the frequency range frex. 
The subplot(211) plot shows the real part of the wavelet in the time domain, 
while the subplot(212) plot shows the power spectrum of the wavelet up to a 
maximum frequency of 80 Hz. The loop pauses after each iteration and waits 
for the user to continue to the next frequency.
%}

figure(5), clf

% hz vector
hz = linspace(0,srate/2,floor(length(wavtime)/2)-1);

for fi=1:nFrex
    
    % create wavelet
    cmw = exp(  2*1i*pi*frex(fi)*wavtime - (wavtime.^2)/(2*s(fi)^2) );
    cmwX = fft(cmw,nConv);
    
    % plot wavelet and its power spectrum
    subplot(211), plot(real(cmw))
    subplot(212), plot(hz,2*abs(cmwX(1:length(hz))));
    
    title([ 'Frequency = ' num2str(frex(fi)) ])
    set(gca,'xlim',[0 80])
    pause % only the user can continue
end

%% try with real EEG data, and illustrate the timetrial

% load EEG data
load ../sampleEEGdata.mat

% create supertrial. What is the size of variable data
% in the next two lines? How could you recover the time-by-trials matrix?
data = squeeze(EEG.data(47,:,:));
data = reshape(data,1,[]);

% create wavelet
wavtime = -2:1/EEG.srate:2;
frex    = 6.5;
s       = 5/(2*pi*frex); % what is the number of cycles here?
csine   = exp(2*1i*pi*frex*wavtime);
gaus    = exp( -(wavtime.^2) / (2*s^2) );
cmw     = csine .* gaus;


% step 1: define convolution parameters
nData = length(data); % or: EEG.pnts*EEG.trials
nKern = length(wavtime);
nConv = nData+nKern-1;
nHfkn = floor(length(wavtime)/2)+1;

% step 2: take FFTs
dataX = fft(data,nConv);
cmwX  = fft(cmw,nConv);

% step 3: normalize kernel
cmwX  = cmwX./max(cmwX); % normalize to get the right units

% step 4: point-wise multiply and take iFFT
as = ifft( dataX.*cmwX );

% step 5: trim wings
as = as(nHfkn:end-nHfkn+1);

% step 5.5: check the size of variable as. Then try plotting it before reshape:
% plot3(1:length(as),real(as),imag(as))

% new step 6: reshape back to time-by-trials
as = reshape(as,EEG.pnts,EEG.trials);


% now plot!
figure(6), clf

subplot(211)
plot(EEG.times,abs(as)*2)
xlabel('Time (ms)'), ylabel('Power')
title([ 'Power at ' num2str(frex) ' Hz from all trials' ])

subplot(212)
plot(EEG.times,mean( abs(as)*2 ,2))
xlabel('Time (ms)'), ylabel('Power')
title([ 'Average power at ' num2str(frex) ' Hz over all trials' ])


% and plot ITPC
figure(7), clf

subplot(211)
plot(EEG.times,angle(as))
xlabel('Time (ms)'), ylabel('Phase angles (rad.)')
title([ 'Phase time series at ' num2str(frex) ' Hz from all trials' ])

subplot(212)
plot(EEG.times,abs( mean(exp(1i*angle(as)),2) ))
xlabel('Time (ms)'), ylabel('Phase clustering')
title([ 'Phase clustering time series at ' num2str(frex) ' Hz over trials' ])

%% edge effects in time-frequency analyses
%{
This code generates a plot showing the time-frequency analysis of a signal using Morlet wavelets. 
The signal is a step function with a duration of 0.8 seconds and a sampling rate of 1000 Hz. 
The Morlet wavelets are created with frequencies ranging from 1 to 100 Hz and a width that 
increases with frequency. The result is a time-frequency plot that shows the power of the 
signal in each frequency band over time.

One issue with time-frequency analysis using wavelets is edge effects, which can lead to 
artifacts at the beginning and end of the signal. In this code, the signal is zero-padded to 
avoid edge effects, but this can introduce additional artifacts if the signal is not zero-mean. 
The code also normalizes the wavelet power spectrum to avoid artifacts due to different 
magnitudes of wavelets at different frequencies.
%}
srate = 1000;
sigtime = 0:1/srate:2;
signal = zeros(size(sigtime));
signal(dsearchn(sigtime',.6):dsearchn(sigtime',1.4)) = 1;

nData = length(sigtime);
nKern = length(wavtime);
nConv = nData+nKern-1;
nHfkn = floor(length(wavtime)/2)+1;

as = ifft( fft(cmw,nConv) .* fft(signal,nConv) );
as = as(nHfkn:end-nHfkn+1);


nFrex = 50;
frex  = linspace(1,srate/10,nFrex);
s     = linspace(4,12,nFrex) ./ (2*pi.*frex);
sigX  = fft(signal,nConv);
tf    = zeros(nFrex,length(signal));

for fi=1:nFrex
    
    % create Morlet wavelet
    cmw = exp(  2*1i*pi*frex(fi)*wavtime - (wavtime.^2)/(2*s(fi)^2) );
    
    % compute its FFT      
    cmwX = fft(cmw,nConv);
    cmwX = cmwX./max(cmwX);
    
    % 'meat' of convolution
    as = ifft( sigX .* cmwX );
    tf(fi,:) = abs( as(nHfkn:end-nHfkn+1) )*2;
end

figure(8), clf

subplot(211)
plot(sigtime,signal)
set(gca,'ylim',[-.05 1.05])

subplot(212)
contourf(sigtime,frex,tf,40,'linecolor','none')
set(gca,'clim',[0 .1])

%% reflection
%{
This code defines a signal g using a Gaussian function, and then generates a 
reflected version of the signal reflex by concatenating the reversed g to the original g, 
and then concatenating that with another reversed g. The resulting reflex signal has the 
same shape as g, but is reflected around its midpoint.

The purpose of generating the reflected signal is to avoid edge effects when performing 
time-frequency analyses, by allowing for a "buffer" of data on either side of the signal. 
By reflecting the signal and then concatenating it with itself, we effectively create a 
"mirror image" of the signal that can be used to "wrap around" the edges of the original signal.

In practical terms, we might use the reflected signal reflex as input to a time-frequency analysis 
function, in place of the original signal g. This would help to mitigate edge effects and ensure 
that the analysis is not unduly influenced by the boundary conditions of the signal.
%}
t = -1:1/srate:3;
g = [ diff( exp(-t.^2/.1) ) 0];
reflex = [ g(end:-1:1) g g(end:-1:1) ];

figure(9),clf

% plot original signal
subplot(211)
plot(g)
set(gca,'xlim',[0 numel(reflex)]-numel(g))

% plot reflected signal
subplot(212)
plot(length(g)+1:2*length(g),g,'ko')
hold on
plot(reflex)
set(gca,'xlim',[0 numel(reflex)])
legend({'original';'reflected'})

%% short-time FFT
%{
This code implements a short-time Fourier transform (STFT) to analyze a chirp signal. 
The signal is created using sine waves with frequencies that change linearly over time (i.e., a chirp). 
The STFT is a way of analyzing how the frequency content of a signal changes over time by breaking 
the signal into smaller time windows and performing a Fourier transform on each window.

The code begins by defining the parameters of the chirp signal and then creating it using the sine function. 
Next, the parameters for the STFT are defined, including the width of the FFT window in milliseconds, 
the number of time steps to analyze, and the indices for each time step. The Hz vector is defined for the FFT window size, 
and a Hann window is created to use during each FFT. The STFT is then computed using a loop over the time steps, 
where the data within the FFT window is extracted, the Hann window is applied, and the FFT is computed. 
Finally, the resulting power spectrum is stored in a matrix for plotting.

The resulting plot shows the original chirp signal in the top subplot and the resulting STFT in the bottom subplot. 
The STFT is represented as a contour plot with time on the x-axis, frequency on the y-axis, and power represented by color. 
The plot shows how the frequency content of the chirp signal changes over time, with the power concentrated in a narrow 
frequency range that gradually increases over time.
%}
% re-create chirp signal
srate = 1000;
sigtime = 0:1/srate:10;
n = length(sigtime);

f  = [2 10];
ff = linspace(f(1),mean(f),n);
signal = sin(2*pi.*ff.*sigtime);



% STFFT parameters
fftWidth_ms = 1000; % FFT width in ms

% convert width from ms to indices
fftWidth = round(fftWidth_ms/(1000/srate)/2);
Ntimesteps = 50; % number of time widths
ct = round(linspace(fftWidth+1,n-fftWidth,Ntimesteps));

% figure out Hz vector for our FFT width
hz = linspace(0,srate/2,fftWidth-1);

hwin = .5*(1-cos(2*pi*(1:fftWidth*2)/ (fftWidth*2-1)));

tf = zeros(length(hz),length(ct));
for ti=1:length(ct)
    tdat = signal(ct(ti)-fftWidth:ct(ti)+fftWidth-1);
    x = fft(hwin.*tdat)/fftWidth;
    tf(:,ti) = 2*abs(x(1:length(hz)));
end


figure(9), clf
subplot(211), plot(sigtime,signal)
set(gca,'ylim',[-1.1 1.1])

subplot(212)
contourf(ct,hz,tf,10,'linecolor','none')
set(gca,'ylim',[0 max(f)*1.5],'clim',[0 1])

%% baseline normalization using dB

time = -2:1/100:5; % 100 Hz sampling rate
basetime = [-1.5 -.5];
baseidx = dsearchn(time',basetime');

mothersig = sin(2*pi*5* (time-3) )./(time-3);
mothersig(~isfinite(mothersig)) = max(mothersig);

figure(10), clf
subplot(311)
plot(time,mothersig)
xlabel('Time (s)'), ylabel('Amplitude')
title('"Mother" signal')


% create two children signals by adding different DC offsets
sig1 = mothersig + 2000;
sig2 = mothersig +  200;

% compute the average per signal
basePow1 = mean(sig1(:,baseidx(1):baseidx(2)),2);
basePow2 = mean(sig2(:,baseidx(1):baseidx(2)),2);

% dB-normalization
sig1DB = 10*log10( bsxfun(@rdivide,sig1,basePow1) );
sig2DB = 10*log10( bsxfun(@rdivide,sig2,basePow2) );

% plot the two signals in the same plot. Easy to compare them?
subplot(312)
plot(time,sig1,time,sig2)
xlabel('Time (sec.)'), ylabel('Raw amplitude')

% and show the normalized signals
subplot(313)
plot(time,sig1DB,time,sig2DB)
legend({'sig1';'sig2'})
xlabel('Time (sec.)'), ylabel('dB amplitude')
title('Normalized signals')

%% plotting time-domain response on top of TF plot

load sampleEEGdata.mat

% extract time-frequency power here...

%% produce a contour plot of TF data
%{
The code below demonstrates baseline normalization using dB. 
It creates two signals with different DC offsets, computes their average baseline power, 
and then normalizes the signals by dividing them by the baseline power and taking the 
logarithm (in dB units) of the result. The normalized signals are then plotted in the bottom 
panel, showing that they have the same baseline level.

dB stands for decibel, which is a unit used to express the ratio of two values, usually the 
ratio of two physical quantities such as power or amplitude. The decibel scale is logarithmic, 
which means that a difference of 1 dB corresponds to a ratio of the two values of 10^(1/10) 
or approximately 1.259. dB is often used in signal processing to express the relative power or 
amplitude of a signal with respect to a reference level. In this context, dB can be used for 
normalization or to compare the amplitudes of two signals.
%}
figure(11), clf
contourf(EEG.times,frex,tf,40,'linecolor','none')
set(gca,'ydir','normal','xlim',[-300 1300],'clim',[-3 3])
xlabel('Time (ms)'), ylabel('Frequency (Hz)')

% now you can rescale the ERP (time-domain trial average)
erp = squeeze(mean(EEG.data(23,:,:),3));
erp = (erp-min(erp))./max(erp-min(erp));
yscale = get(gca,'ylim');
erp = erp*(yscale(2)-yscale(1))+yscale(1);

hold on
plot(EEG.times,erp,'k','linew',2)

%% end.
