%% Matlab for Brain and Cognitive Scientists
% Matlab code for Chapter 20
% Mike X Cohen
% 
% This code accompanies the book 
%      "Matlab for Brain and Cognitive Scientists" (MIT Press, 2017).
% Using the code without following the book may lead to confusion, 
% incorrect data analyses, and misinterpretations of results. 
% Mike X Cohen assumes no responsibility for incorrect use of this code. 

%% running-mean filter
%{
This code implements a running-mean filter on a simulated time-series signal. 
The signal is generated as a sinusoidal wave with a frequency of 8 Hz and a random noise component. 
The running-mean filter is then applied to the signal using a filter window of 7 samples.
%}
srate = 1000;
time  = 0:1/srate:6;
ampl  = interp1(0:5,rand(6,1),time,'spline');
%{
interp1 is a function in MATLAB and Octave that performs one-dimensional interpolation of data. 
In the provided code, the interp1 function is used to interpolate a random amplitude signal ampl from 
a given set of time points (0:5) to a new set of time points (time).

The first argument of interp1 is the original set of times for which we have the ampl values (in this case, 0:5). 
The second argument is the original set of ampl values (a random vector of size 6x1). The third argument specifies 
the new set of times for which we want to interpolate the ampl values (in this case, time variable). 
Finally, the string argument 'spline' specifies the type of interpolation to use (in this case, cubic spline interpolation).
%}
f     = 8; % Hz because time is in seconds
noise = .3*randn(size(time));
signal = ampl .* sin(2*pi*f*time) + noise;
%{
In this line, a signal is generated by multiplying a sinusoidal wave with a random amplitude function ampl and 
adding some noise.
%}
k = 7;
%{
In the given code, k is used as the window size for the running-mean filter. Specifically, k is the number of 
samples on each side of the current sample that are included in the average calculation. For example, 
if k is 7, then the running-mean filter will use 15 samples in total to calculate the mean at each time 
point (7 samples to the left, the current sample, and 7 samples to the right).
%}
filtsig = zeros(size(signal));
%{
In this line, a zero array of the same size as the original signal is created. This array will be used to store the filtered signal.
%}
for i=k+1:length(time)-k
    filtsig(i) = mean(signal(i-k:i+k) );
end
%{
This is a for loop that iterates over the time points in the signal from k+1 to length(time)-k. 
The reason why the loop starts at k+1 is that it's not possible to compute the mean of k 
samples before the first time point, and it's not possible to compute the mean of 
k samples after the last time point.

For each iteration of the loop, the mean of the signal values at the current time point 
and the k time points before and after the current time point are computed using the 
mean function, and the result is stored in the filtsig variable at the current time point index i.
%}

figure(1), clf
plot(time, signal, time, filtsig)
legend({'original';'filtered'})
title('Illustration of running-mean filter')

%% running-median filter
%{
This code demonstrates the implementation of a running-median filter to remove noise from a noisy signal.

The difference between the running-mean and running-median codes is that the running-mean code computes 
the mean of the signal over a certain number of adjacent samples (defined by the k parameter), 
while the running-median code computes the median of the signal over the same window.

In the running-median code, the sort function is used to sort the samples within the window, 
and the median is then defined as the middle sample if the number of samples in the window is odd, 
or the average of the two middle samples if the number of samples is even.

Another difference between the two codes is that the running-median code adds 
a noise signal to the input signal (signal variable) to better illustrate the effects of the filtering.
%}
srate = 1000;
time  = 0:1/srate:6;
ampl  = interp1(0:5,rand(6,1),time,'spline');
f     = 8; % Hz because time is in seconds
noise = 1000*isprime(1:length(time));
%{
In the code, isprime(1:length(time)) generates a logical vector of the same size as time, with elements 
true or false indicating whether each element of 1:length(time) is a prime number or not.

Multiplying this logical vector by 1000 results in a vector of the same size with values of 1000 for prime 
numbers and 0 for non-prime numbers. This creates a "noisy" signal with values that randomly spike to 
1000 at the locations of prime numbers in the time vector.
%}
signal = ampl .* sin(2*pi*f*time) + noise;

k = 7;
[filtMean, filtMedian] = deal( zeros(size(signal)) );
%{
This line initializes two variables filtMean and filtMedian as zero vectors of the same size as the signal vector, 
using the zeros function. These variables will be used to store the filtered signals, one filtered using a 
running-mean filter and the other using a running-median filter.

Here, the deal function is used to assign the same zero vector to both variables simultaneously.
%}
for i=k+1:length(time)-k
    filtMean(i) = mean( signal(i-k:i+k) );
    
    sortnums = sort( signal(i-k:i+k) );
    filtMedian(i) = sortnums(k);
end
%{
For the median, the sort() function is used to sort the 2*k+1 points in ascending order. 
The kth point in this sorted array is selected as the median, and is stored in the corresponding element of filtMedian.

Sorting the window of data points is necessary to identify the middle value(s) in the window, which is the median.
%}

figure(2), clf

subplot(221)
plot(time,signal)
set(gca,'xlim',[1 2],'ylim',[-1.5 1.5])
title('Signal')
xlabel('Time (s)'), ylabel('Amplitude')

subplot(222)
plot(time,filtMean)
title('Mean-filtered')
xlabel('Time (s)'), ylabel('Amplitude')
set(gca,'xlim',[1 2],'ylim',[-50 350])


subplot(223)
plot(time,filtMedian)
title('Median-filtered')
xlabel('Time (s)'), ylabel('Amplitude')
set(gca,'xlim',[1 2],'ylim',[-1.5 1.5])

zoom on

%% edge in the frequency domain and its time-domain representation
%{
The code is plotting the frequency-domain and time-domain representation of an edge or a sharp point. 
In the frequency domain, the code sets the value of the 5th percentile of the data in the vector X to 1, 
and all other elements are 0. This creates a sharp point in the frequency domain. 
The code then uses the ifft function to perform the Inverse Fast Fourier Transform
 to convert the frequency-domain signal back to the time-domain signal. 
Finally, the plot function is used to plot the time-domain signal, and the set 
function is used to set the y-axis limits to [-1.1 1.1]. The signal in the time-domain is 
normalized by the length of X using the expression real(ifft(X))*length(X).
%}
% sharp point (a.k.a. impulse response function, a.k.a. delta function)
N = 400;
X = zeros(N,1);
X(round(N*.05)) = 1;        % This line of code sets the value of the 5th percentile of the data in the vector X to 1.

figure(3), clf
subplot(221)
plot(X)
set(gca,'ylim',[-.01 1.05])
xlabel('Frequency (a.u.)'), ylabel('amplitude')

subplot(223)
plot(real(ifft(X))*length(X))
%{
X represents the frequency-domain signal, and ifft is the Inverse Fast Fourier Transform function, 
which is used to convert the frequency-domain signal back to the time-domain signal. 
The real function is used to extract the real part of the signal. Finally, the signal is scaled 
by the length of X to normalize the amplitude.
%}
set(gca,'ylim',[-1.1 1.1])
xlabel('Time (a.u.)'), ylabel('amplitude')


% box-car
X = zeros(N,1);
X(10:30) = 1;

subplot(222)
plot(X)
set(gca,'ylim',[-.01 1.05])
xlabel('Frequency (a.u.)'), ylabel('amplitude')

subplot(224)
plot(real(ifft(X))*length(X))
set(gca,'ylim',[-15 15])
xlabel('Time (a.u.)'), ylabel('amplitude')


%% empirical FWHM of a Gaussian
%{
You want the filter to be narrow to maximize frequency specificity, but you want the filter to have
gentle slopes rather than sharp edges, although this reduces frequency specificity.
What shape might come into mind when you think about a narrow
peak that has smoothly decaying sides?

If you weren�ft already thinking that the answer is a Gaussian, then I�fm
sure you are thinking it now. Perhaps you remember from the previous
chapter that Morlet wavelets have a Gaussian shape in the frequency
domain. Here we will construct a similar filter, except we define the kernel
in the frequency domain instead of in the time domain.
%}
% create a gaussian
x = -4:.1:4;
gaus = exp(-(x.^2));
gaus = gaus./max(gaus);
%{
FWHM is the distance between the 50% amplitude points before to after
the peak (figure 20.4). If the Gaussian is normalized to have an amplitude
of 1.0 (which can be obtained by dividing by the maximum), then the task
is simplified by finding the point before the peak and after the peak closest
to 0.5. We start by creating a Gaussian and maximum-normalizing it.

Note: gaus is not a MATLAB function!
%}

% find index of peak
[~,pidx] = max(gaus);
prepeak  = dsearchn(gaus(1:pidx)',.5);
postpeak = pidx-1+dsearchn(gaus(pidx:end)',.5);
%{
Our next step is to find the indices before and after the peak that are
closest to 0.5.

gaus(1:pidx)' returns the values of the gaus vector from the beginning 
to the index of the peak, while gaus(pidx:end)' returns the values of the 
gaus vector from the index of the peak to the end. These values are used 
as input to dsearchn to find the indices of the values closest to 0.5. pidx-1 is 
added to dsearchn(gaus(pidx:end)',.5) to account for the fact that the index 
of the peak in the original gaus vector is not 1.
%}

% plot gaussian
figure(4), clf
plot(x,gaus,'k-o')
hold on

% plot the empirical points closest to 50% amplitude
plot(x(prepeak),gaus(prepeak),'ro')
plot(x(postpeak),gaus(postpeak),'ro')
plot([ x(prepeak) x(postpeak)],[gaus(prepeak) gaus(postpeak)],'r--')
plot(get(gca,'xlim'),[.5 .5],'k:')
ylabel('Normalized amplitude')

% so what is the FWHM in Hz?
%{
The Gaussian narrow-band filter has three advantages. It is easy to construct
and implement in MATLAB; it produces minimal ripple artifacts in
the time-domain reconstruction, because there are no edges in the frequency
domain; and it is a non-causal filter, meaning it does not produce
any phase distortions in the result (although it does mean that there is
temporal leakage �gbackward�h in time as well as forward in time). There are
two disadvantages. The result is real-valued, so you�fll need to apply the
Hilbert transform if you want to extract power and phase (more about this
in a later section). Perhaps the main disadvantage is that you cannot control
the shape of the filter. The Gaussian filter is?I feel silly writing this but
it segues to the next section?shaped like a Gaussian. If you want one filter
that spans, say, 30?60 Hz, a Gaussian filter is a suboptimal choice.
%}
%% create a frequency-domain Gaussian with specified FWHM in hz

fwhm     =  5; % in Hz
centfreq = 14; % center frequency, also in Hz


% set parameters
srate = 1000;
N = 4000;

hz = linspace(0,srate,N);  % frequencies
s  = fwhm*(2*pi-1)/(4*pi); % normalized width 
%{
The line s = fwhm*(2*pi-1)/(4*pi); computes the standard deviation of the Gaussian curve in radians, 
where fwhm is the full-width at half maximum of the curve.
%}
x  = hz-centfreq;          % shifted frequencies 
gx = exp(-.5*(x/s).^2);    % gaussian 
gx = gx./max(gx);          % gain-normalized 



% compute empirical frequency and standard deviation
idx = dsearchn(hz',centfreq);
%{
dsearchn is a function in MATLAB that finds the indices of the closest matching values in an array. 
Here, it is being used to find the index of the frequency value closest to the center frequency 
specified in the variable centfreq. The apostrophe in hz' is a MATLAB notation for taking the 
transpose of the hz vector to ensure that it is a column vector for use in the dsearchn function.
%}
empVals(1) = hz(idx);
%{
This line of code finds the index of the frequency closest to the center frequency centfreq and 
assigns the corresponding frequency value to empVals(1).
%}
% find values closest to .5 after MINUS before the peak
empVals(2) = hz(idx-1+dsearchn(gx(idx:end)',.5)) - hz(dsearchn(gx(1:idx)',.5));
%{
This line of code computes the empirical FWHM in Hz. It first finds the index of the peak of the 
Gaussian using the previously computed idx variable. Then it finds the indices of the values in gx 
closest to 50% of the peak amplitude, separately for the values to the left and right of the peak. 
It subtracts the frequency at the index closest to 50% on the left from the frequency at the index 
closest to 50% on the right to get the empirical FWHM in Hz.
%}

% plot
figure(5), clf
plot(hz,gx,'o-')
hold on
plot([hz(dsearchn(gx(1:idx)',.5)) hz(idx-1+dsearchn(gx(idx:end)',.5))], [gx(dsearchn(gx(1:idx)',.5)) gx(idx-1+dsearchn(gx(idx:end)',.5))],'k--')
%{
This line of code plots a dashed line connecting the two points where the Gaussian curve intersects the half-maximum level. 
It uses the dsearchn function to find the indices of the frequency values corresponding to these points. 
The dsearchn function finds the index of the nearest element in the second input argument (0.5 in this case) 
within the vector or matrix specified in the first input argument (gx(1:idx) and gx(idx:end) in this case). 
The indices are then used to retrieve the corresponding frequency values from the hz vector and plot the 
dashed line using plot.
%}
set(gca,'xlim',[max(centfreq-fwhm*2,0) centfreq+fwhm*2],'ylim',[-.01 1.01]);

title([ 'Requested: ' num2str(centfreq) ', ' num2str(fwhm) ' Hz; Empirical: ' num2str(empVals(1)) ', ' num2str(empVals(2)) ' Hz' ])
xlabel('Frequency (Hz)'), ylabel('Amplitude gain')

%% now we apply the filter to a chirp

% re-create chirp signal
% same length as Gaussian created in previous cell!
sigtime = 0:1/srate:4-1/srate;
n = length(sigtime);

f  = [2 10];
ff = linspace(f(1),mean(f),n);
signal = sin(2*pi.*ff.*sigtime);
%{
The chirp signal is defined in the above code block.

A chirp signal is a signal whose frequency varies with time. 
The frequency can increase, decrease or sweep through a range of 
frequencies over time. This type of signal is commonly used in radar, 
sonar, and in communication systems. A common example of a chirp 
signal is the sound of a bird's chirp, which increases or decreases in pitch over time.
%}

fwhm     =  3; % in Hz
centfreq =  8;


N = length(sigtime);

hz = linspace(0,srate,N);  % frequencies
s  = fwhm*(2*pi-1)/(4*pi); % normalized width 
x  = hz-centfreq;          % shifted frequencies 
gx = exp(-.5*(x/s).^2);    % gaussian 
gx = gx./max(gx);          % gain-normalized 


% apply the filter
filtsig = 2*real( ifft( fft(signal).*gx ) );

% and plot
figure(6), clf
h=plot(sigtime,signal);
hold on
plot(sigtime,filtsig,'k','linew',2)

%% setting up an FIR filter
%{
In this code block, we are setting up a Finite Impulse Response (FIR) filter. 
The filter's frequency response is defined by its shape and frequency values. 
Initially, the shape and frequency values are not well-defined and the frequency 
response plot looks pretty bad. We then define the transition zone and frequency 
boundaries, and update the shape and frequency values to create a better frequency response plot.
%}
% sampling rate
srate = 500;

% shape of the frequency response
fshape = [ 0 1 1 0 ];
frex = [ 0 12 18 srate/2 ]./srate/2;

% plot it.
figure(7), clf

subplot(211)
plot(frex,fshape,'k-o','linew',3,'markersize',8)
xlabel('Frequency (fraction of Nyquist)'), ylabel('Gain')
set(gca,'ylim',[-.02 1.02])

% That doesn't look very good. Actually, it looks pretty awful.
% Let's add some more parameters...

tz = .15; % transition zone, in percent
fbnd = [ 12 18 ]; % frequency boundaries
fshape = [ 0 0 1 1 0 0 ]; % note the difference between this and the previous shape
frex = [ 0 fbnd(1)*(1-tz) fbnd fbnd(2)*(1+tz) srate/2 ];
frex = frex./(srate/2); % norm. to Nyquist

% plot the improved filter response
subplot(212)
plot(frex,fshape,'k-o','linew',3,'markersize',8)
xlabel('Frequency (Hz)'), ylabel('Gain')
set(gca,'ylim',[-.02 1.02])

%% apply filter to linear chirp
%{
This code applies a finite impulse response (FIR) filter to a linear chirp signal. 
The chirp signal is created using Matlab's chirp function, and the FIR filter is 
created using Matlab's firls function. The filter is designed to have a band-pass shape, 
allowing only frequencies within a certain range to pass through while attenuating frequencies 
outside that range. The filtered signal is obtained using the filtfilt function, which applies the filter 
forward and backward to minimize phase distortion. The code then plots the filter kernel in the time
and frequency domains, as well as the original and filtered signals.

A linear chirp is a signal whose frequency changes linearly with time. It starts at a certain frequency 
and increases (or decreases) linearly with time. A linear chirp is a type of signal that is commonly used in 
signal processing and communication applications. It is often used in radar, sonar, and ultrasound imaging 
to analyze the properties of objects or materials.
%}

% this time we'll use Matlab's chirp function
srate  = 500;
time   = 0:1/srate:15;
frange = [5 20];
x  = chirp(time, frange(1), time(end), frange(2));  % for octave: 'chirp' is in the signal package
%{
chirp is a function in MATLAB that generates a linearly swept frequency signal. 
The first input argument time specifies the time axis for the signal, frange(1) and frange(2) 
specify the starting and ending frequencies of the sweep, respectively. The resulting signal x is a 
time-varying sinusoid whose frequency increases linearly from frange(1) to frange(2) over the 
time interval specified by time. The chirp signal is often used in signal processing and communication 
applications, such as radar and sonar.
%}

% setup the FIR filter again
tz     = .15; % transition zone, in percent
fbnd   = [ 13 15 ]; % freq boundaries
fshape = [ 0 0 1 1 0 0 ];
frex   = [ 0 fbnd(1)*(1-tz) fbnd fbnd(2)*(1+tz) srate/2 ];
frex   = frex./(srate/2); % norm. to Nyquist
%{
These lines of code define the filter specifications for a Finite Impulse Response (FIR) filter. 
The filter is designed to have a passband between 13 Hz and 15 Hz, which is defined by the fbnd variable.

The fshape variable specifies the desired shape of the frequency response of the filter. 
In this case, it is a rectangle with a height of 1 in the passband and 0 in the stopbands.

The frex variable defines the frequency values that correspond to the shape specified by fshape. 
These values are normalized to the Nyquist frequency (i.e., half the sampling rate).

The tz variable defines the transition zone width as a percentage of the passband width. 
The transition zone is the region between the passband and stopbands where the filter 
response gradually transitions from the passband to the stopbands. The frex variable is 
then modified to include the transition zone by extending the passband boundaries.
%}

% filter order
orderfact = 3;
filtord = round( (orderfact*1000/fbnd(1)) / (srate/1000) );
%{
filtord is the filter order, which determines the length of the filter kernel used to filter the data. 
It is calculated based on the desired frequency band of the filter (fbnd) and the desired order factor (orderfact).
%}

% create the filter kernel using Matlab's firls function
filtkernel = firls(filtord,frex,fshape); % for octave: also in signal package

% filter the data
fdata = filtfilt(filtkernel,1,x); % for octave: also in signal package
%{
firls() is a MATLAB function that is used to design a linear-phase FIR filter using least squares approximation. 
The filter coefficients are obtained by fitting a piecewise constant function to the desired frequency response 
 by fshape over the frequency bands specified by frex, with the transition band width given by tz. 
The filter order is determined by filtord.

In this particular code, firls() is used to design a band-pass filter with frequency response fshape in the frequency 
band specified by fbnd with a transition zone tz of 15%. The normalized frequencies of the filter are stored in frex, 
which is normalized to the Nyquist frequency. The filter order is determined by filtord, which is calculated based on 
orderfact, a factor used to determine the filter order.
%}

figure(8), clf
subplot(221)
kernelTime = linspace(1,orderfact*1000/fbnd(1),filtord+1);
plot(kernelTime, filtkernel)
set(gca,'xlim', kernelTime([1 end]))
xlabel('Time (ms)')
title('Filter kernel in the time domain')

subplot(222)
hz = linspace(0, srate/2, floor((filtord + 1)/2)+1);
powr = abs(fft(filtkernel)).^2;
plot(hz,powr(1:length(hz)))
hold on
plot(frex*srate/2,fshape,'ro-')
set(gca,'xlim',[0 fbnd(end)*3],'ylim',[-.05 1.1])
legend({'Empirical filter';'''Ideal'' filter'})
xlabel('Frequency (Hz)')
title('Filter kernel in the frequency domain')


subplot(212)
plot(time,x,'m'), hold on
plot(time,fdata,'k','linew',2)
% set(gca,'xlim',[6 9])
pan on

%% hilbert transform

% a random signal and its Fourier coefficients
n = 20;
r = randn(n,1);
rx = fft(r);

% a copy that is multiplied by the complex operator
rxi = rx*1i; 

% find indices of positive and negative frequencies
posF = 2:floor(n/2)+mod(n,2); 
negF = ceil(n/2)+1+~mod(n,2):n; 

% rotate Fourier coefficients 
rx(posF) = rx(posF) + -1i*rxi(posF); 
rx(negF) = rx(negF) +  1i*rxi(negF); 

% take inverse FFT 
hilbert_r = ifft(rx);

figure(9), clf

subplot(311)
plot(r,'o-','markersize',9)
hold on
plot(real(hilbert_r),'k*')
plot(real( hilbert(r) ),'gd','markersize',20) % for octave: hilbert function is in the signal package
legend({'original';'real part of Hilbert';'real part of Matlab hilbert'})
title('Real part')


subplot(312)
plot(abs(hilbert_r).^2,'k*-')
hold on
plot(abs(hilbert(r)).^2,'gd','markersize',20)
xlabel('Frequency (a.u.)'), ylabel('Amplitude (a.u.)')
legend({'power from Hilbert';'power from Matlab hilbert'})
title('Power')


subplot(313)
plot(angle(hilbert_r),'k*-')
hold on
plot(angle(hilbert(r)),'gd','markersize',20)
xlabel('Frequency (a.u.)'), ylabel('Phase (rad.)')
legend({'phases from Hilbert';'phases from Matlab hilbert'})
title('Phases')

%% end.
