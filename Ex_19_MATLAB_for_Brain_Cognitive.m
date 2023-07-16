%% 19.7.1
%{
In the code for figure 19.5, uncenter the Morlet wavelet by defining
wavtime from -2 to +3 seconds (also try other ranges like -2 to +2.1).
How does this change the time-frequency plot, why does this happen,
and what is the important lesson about constructing wavelets?

It moves the time-frequency plot towards the side, depending on the time
range of wavtime variable.

The reason why the contour plot moves to the left when wavtime parameters 
are changed is due to the shift property of the Fourier transform. When the time 
window of the Morlet wavelet is shifted, the frequency content of the wavelet changes 
as well. This is because the wavelet's Fourier transform has a phase factor that is 
dependent on the time parameter. Therefore, when you change the time window of the 
wavelet, you are effectively changing the phase of its Fourier transform. This phase shift 
causes the frequency content of the wavelet to shift in the frequency domain, which in turn 
causes the time-frequency representation to shift in the time domain.

The important lesson here is to match the time parameter of the frequency
of the Morlet waves to that of the signal. 
%}

% Uncentered Morlet Wavelet:
wavtime2 = -2:1/srate:3;
frex = 6.5;
s = 15/(2*pi*frex);
csine = exp(2*1i*pi*frex*wavtime2);
gaus = exp( -(wavtime2.^2) / (2*s^2) );
cmw = csine .* gaus;

sigtime = 0:1/srate:2;
signal = zeros(size(sigtime));
signal(dsearchn(sigtime',.6):dsearchn(sigtime',1.4)) = 1;

nData = length(sigtime);
nKern = length(wavtime2);
nConv = nData+nKern-1;
nHfkn = floor(length(wavtime2)/2)+1;

as = ifft( fft(cmw, nConv) .* fft(signal, nConv) );
as = as(nHfkn:end-nHfkn+1);


nFrex = 50;
frex  = linspace(1,srate/10,nFrex);
s     = linspace(4,12,nFrex) ./ (2*pi.*frex);
sigX  = fft(signal,nConv);
tf    = zeros(nFrex,length(signal));

for fi=1:nFrex
    
    % create Morlet wavelet
    cmw = exp(  2*1i*pi*frex(fi)*wavtime2 - (wavtime2.^2)/(2*s(fi)^2) );
    
    % compute its FFT      
    cmwX = fft(cmw,nConv);
    cmwX = cmwX./max(cmwX);
    
    % 'meat' of convolution
    as = ifft( sigX .* cmwX );
    tf(fi,:) = abs( as(nHfkn:end-nHfkn+1) )*2;
end

figure(1), clf

subplot(211)
plot(sigtime,signal)
set(gca,'ylim',[-.05 1.05])

subplot(212)
contourf(sigtime, frex, tf, 40,'linecolor','none')
set(gca,'clim',[0 .1])

%% 19.7.4
%{
Adapt the code from figure 19.5 by replacing the boxcar function
with a sine wave at 30 Hz. Confirm that the time-frequency result is an
accurate representation. Then decibel-normalize the result. First use a
baseline period of 0 to 0.2 seconds, then use a baseline of 0.2 to 0.5
seconds, and finally use a baseline of 0.8 to 1.2 seconds. Why do the
results look so different for these different baseline periods?

A boxcar function is a type of function that has a constant value between 
two specified limits and is zero outside those limits. It is also known as a 
rectangular function or a pulse function. The boxcar function is commonly 
used in signal processing and digital signal analysis as a windowing function 
to isolate a specific portion of a signal for analysis. The boxcar function has a 
wide bandwidth in the frequency domain, which can cause spectral leakage 
in the FFT of a signal, so other windowing functions such as Hamming or 
Hanning windows are often used instead.

The results look different because the sine wave frequencies change
intervals, which creates changes in the Morlet at the frequency-time domain
as well. 
%}
srate = 1000;
wavtime = -2:1/srate:2;
frex = 6.5;
s = 5/(2*pi*frex);
csine = exp(2*1i*pi*frex*wavtime);
gaus = exp( -(wavtime.^2) / (2*s^2) );
cmw = csine .* gaus;

% Replace boxcar function:
%{
sigtime = 0:1/srate:2;
signal = zeros(size(sigtime));
signal(dsearchn(sigtime',.6):dsearchn(sigtime',1.4)) = 1;
%}
% With a sine wave function:
sigtime = 0:1/srate:0.2;
freq = 30;
signal = sin(2*pi*freq*sigtime);
%Decible normalize the result:
signal = log(10)*(signal);
signal = signal./max(signal(:));

nData = length(sigtime);
nKern = length(wavtime);
nConv = nData+nKern-1;
nHfkn = floor(length(wavtime)/2)+1;

nFrex = 50;
frex  = linspace(1,srate/10,nFrex);
s     = linspace(4,12,nFrex) ./ (2*pi.*frex);
sigX  = fft(signal,nConv);
tf    = zeros(nFrex, length(signal));

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

figure(2), clf

subplot(211)
plot(sigtime,signal)
set(gca,'ylim',[-.05 1.05])

subplot(212)
contourf(sigtime,frex,tf,40,'linecolor','none')
set(gca,'clim',[0 .1])

%% 19.7.9
%{
The power spectrum of the ÅgtriangleÅh signal from figure 18.6 is nearly
uninterpretable. Apply a time-frequency analysis of this signal. Does
this better capture the dynamics in the data?
%}

% define some signal parameters
srate = 1000;
t = 0:1/srate:10;
n = length(t);
f = 3; % frequency in Hz
scales = 3:30; % scales corresponding to each frequency

% create time-increasing amplitude, and stationary amplitude
ampl1 = linspace(1,10,n);
% the next line creates a randomly varying amplitude, FYI
% ampl1 = abs(interp1(linspace(t(1),t(end),10),10*rand(1,10),t,'spline'));
ampl2 = mean(ampl1);

% Create Morlet wavelets
wavelets = cmorwavf(0, 30, scales);

% Create signal
signal1 = ampl1 .* sin(2*pi*f*t);
signal2 = ampl2 .* sin(2*pi*f*t);

% Compute continuous wavelet transform
cwt_result1 = cwt(signal1, wavelets);
cwt_result2 = cwt(signal2, wavelets);
% Compute power
power1 = abs(cwt_result1).^2;
power2 = abs(cwt_result2).^2;
% Create time-frequency plot
figure(6), clf
subplot(121)
h(1) = imagesc(t, 1./scales, power1);
set(gca, 'ydir', 'normal');
hold on
xlabel('Time'), ylabel('Frequency')

subplot(122)
h(2) = imagesc(t, 1./scales, power2);
set(gca, 'ydir', 'normal');
hold on
xlabel('Time'), ylabel('Frequency')

%% 19.7.11
%{
There is a MATLAB function called chirp, which creates a chirp.
On the basis of the help information, create a chirp using the same
parameters as for the chirp we created manually in this chapter. Plot
them against each other. Do they match? (The answer is no.) How
much are they off, and how can you adjust your code to match the
output of chirp?
%}