%% Matlab for Brain and Cognitive Scientists
% Matlab code for Chapter 21
% Mike X Cohen
% 
% This code accompanies the book 
%      "Matlab for Brain and Cognitive Scientists" (MIT Press, 2017).
% Using the code without following the book may lead to confusion, 
% incorrect data analyses, and misinterpretations of results. 
% Mike X Cohen assumes no responsibility for incorrect use of this code. 

%% basic RMS on empirical data
% This is an example of basic RMS on empirical data.

load(['C:\Users\dun81\OneDrive\Documents\Dun II\ET Online\Library Online\Psychology\Computational Neuroscience\MATLAB for Brain and Cognitive Scientists - The MIT Press (2017)\MATLAB_BrainAndCognitiveScientists-main\ch21\eeg1.mat'])

% compute RMS. pretty simple.
%{
The root mean square (RMS) is computed for each trial. 
This is done by first squaring the values of each trial (eeg.^2), 
computing the mean of each trial (mean( eeg.^2 ,1)), and then taking 
the square root of the means (sqrt( mean( eeg.^2 ,1) )). 
The result is a vector rmsx containing the RMS for each trial.
%}
rmsx = sqrt( mean( eeg.^2 ,1) );

% plot the time-domain average response
%{
The time-domain average response is a way to visualize the average response of a signal across time. 
It is obtained by averaging the signal across trials or repetitions. In the case of EEG data, for example, 
the time-domain average response shows the average electrical activity recorded at different electrodes across time. 
This can reveal important features of the signal, such as the timing and amplitude of peaks and troughs, and can be used 
to identify event-related potentials (ERPs) associated with specific cognitive processes or events.
%}
figure(1), clf
subplot(211)
plot(EEGtime, mean(eeg, 2))
xlabel('Time (ms)')
ylabel('Voltage (\muV)')
set(gca,'xlim', EEGtime([1 end]))

% bar plot of the RMS per trial. Notice anything striking?
subplot(212)
bar(rmsx, 'histc')
ylabel('RMS')
xlabel('Trials')
set(gca,'xlim',[0 size(eeg,2)+2])
%{
To make things interesting, IÅfve added noise to one trial. The noise is randomly distributed
around zero, making it almost unnoticeable when inspecting the trial
average. But the RMS plot clearly reveals the outlier trial. For this reason,
RMS can facilitate the initial stages of data cleaning and quality control.
%}

%% RMS over space per time point
%{
This code is trying to demonstrate how to compute the root mean square
(RMS) over space per time point of EEG data. It shows how to compute the
RMS for each trial and the RMS for the average of all trials, and how to
plot both values over time. The purpose is to illustrate how the RMS can be
used to summarize the amplitude of the EEG signal across space and time.
%}
load(['C:\Users\dun81\OneDrive\Documents\Dun II\ET Online\Library Online\Psychology\Computational Neuroscience\MATLAB for Brain and Cognitive Scientists - The MIT Press (2017)\MATLAB_BrainAndCognitiveScientists-main\ch21\eeg1.mat'])

% compute RMS per trial or on trial-average
rmsx_trials = squeeze( sqrt( mean( eeg.^2, 1) ) );
rmsx_averag = squeeze( sqrt( mean( mean(eeg, 3).^2, 1) ) );
%{
In this code, eeg is a 3D matrix of EEG data, where the first dimension corresponds to time, 
the second dimension corresponds to electrodes, and the third dimension corresponds to trials. 
The code computes the RMS (root mean square) value for each electrode and trial separately.

Specifically, mean(eeg.^2, 1) computes the mean squared value for each electrode and trial, 
resulting in a 2D matrix of size num_electrodes x num_trials. The squeeze function is then used 
to remove the singleton dimension (i.e., the first dimension), resulting in a 1D array of size num_trials. 
Finally, sqrt is applied element-wise to obtain the RMS value for each trial.

Similarly, mean(mean(eeg, 3).^2, 1) computes the mean squared value across trials for each electrode, 
resulting in a 1D array of size num_electrodes. Again, squeeze is used to remove the singleton dimension, 
and sqrt is applied element-wise to obtain the RMS value for each electrode.
%}

% and plot
figure(2), clf
plot(EEGtime, mean(rmsx_trials, 2)); 
hold on; 
plot(EEGtime(1:length(rmsx_averag)), rmsx_averag);          % I had to re-write this part b/c it wasn't working due to the plots having different sizes. I changed the lengths of the first variable to be equal to the second variable.
xlim([-1000, max(EEGtime(1:length(rmsx_averag)))]);         % Because the length of the x-variable has been changed, the x-axis limit needs to be changed to reflect this as well. 
legend({'All trials'; 'Trial average'})
xlabel('Time (ms)'), ylabel('\muV^2')

%% DFA on human task data

% load data
load(['C:\Users\dun81\OneDrive\Documents\Dun II\ET Online\Library Online\Psychology\Computational Neuroscience\MATLAB for Brain and Cognitive Scientists - The MIT Press (2017)\MATLAB_BrainAndCognitiveScientists-main\ch21\behaviorData.mat'])

N = length(perfdat);

% setup parameters
nScales = 20;
ranges  = [1 400]*srate;
scales  = ceil(logspace(log10(ranges(1)), log10(ranges(2)), nScales));
rmses   = zeros(size(scales));
%{
In this code, the variables nScales, ranges, scales, and rmses are being defined to set up the parameters 
for the Detrended Fluctuation Analysis (DFA) computation on the perfdat variable.

nScales is the number of scales that will be used to compute the DFA, and ranges is a two-element vector 
containing the minimum and maximum scale values. scales is defined by computing a logarithmically-spaced 
set of nScales values between the minimum and maximum scale values, using the logspace function. 
Finally, rmses is initialized as a vector of zeros with the same size as scales, which will be used to store the 
computed RMS values at each scale.
%}

% plot the data
figure(3), clf
subplot(211)
plot(timevec/60, perfdat)
xlabel('Time (minutes)'), ylabel('Performance')
set(gca,'ytick',[-1 1],'yticklabel',{'error';'correct'},'ylim',[-1.1 1.1],'xlim', timevec([1 end])/60)
%{
This code generates a figure with two subplots. The top subplot shows the performance data over time in minutes, 
with the y-axis labeled as "error" or "correct" depending on whether the participant made an error or a correct response. 
The x-axis represents time in minutes. The bottom subplot is currently empty.

perfdat is a variable that contains performance data. It is loaded from the behaviorData.mat file and is used in the code to demonstrate 
self-similar behavior and perform DFA analysis on it. 

timevec is a vector containing time points at which the performance data in perfdat was sampled. It is likely that timevec has the same 
number of elements as perfdat, and each element corresponds to a specific time point.
%}

% integrate data
perfdat = cumsum(perfdat(:)-mean(perfdat));

% and show that for comparison
subplot(212)
plot(timevec/60, perfdat)
set(gca,'xlim', timevec([1 end])/60)
xlabel('Time (minutes)')


% qualitative demonstration of self-similar behavior
figure(4), clf
for i=1:4
    subplot(4,1, i)
    plot(timevec/60, detrend(perfdat))
end
%{
detrend is a MATLAB function that removes a linear trend from a data vector or matrix along a specified dimension. 
It calculates the least-squares fit of a straight line to the data and subtracts it. This can be useful in removing systematic 
trends or biases in data, leaving behind only the fluctuations that are of interest.
%}
subplot(411), set(gca,'xlim',[0 40])
subplot(412), set(gca,'xlim',[25 29])
subplot(413), set(gca,'xlim',[27 27.4])
subplot(414), set(gca,'xlim',[27.2 27.24])

%% an aside on detrend...
%{
This code generates a random noise signal with a linear trend and then applies the detrend() function in MATLAB 
to remove the linear trend from the signal. It then plots both the original and detrended signals for comparison.
%}
% 200 points, because why not?
n = 200;

% random noise plus a linear trend
d = randn(1,n) + linspace(-3,10,n);

figure(5), clf
plot(1:n,d, 1:n, detrend(d) )
legend({'Original';'Detrended'})

%% compute RMS over different time scales
%{
This code segment computes the root mean square (RMS) over different time scales (or windows) of the input data perfdat. 
The scales variable determines the size of each window by defining the range of values that is log-spaced (logspace) and then 
the ceil function is applied to ensure integer values. The code then loops through the range of scales and calculates the RMS for 
each window using the detrend function to remove linear trends within each window. Finally, the code fits a linear model to the 
log-transformed scale and RMS data to quantify the scaling exponent using a method called detrended fluctuation analysis (DFA).
%}

for scalei = 1:length(scales)
    
    % epoch
    n = floor(N/scales(scalei)); % number of epochs
    epochs = reshape( perfdat(1:n*scales(scalei)) ,scales(scalei),n)';
    
    % detrend
    depochs = detrend(epochs')';
    % note the effect of detrending:
    %clf, subplot(211), plot(epochs'), subplot(212), plot(depochs')
    
    % RMS and average
    rmses(scalei) = mean( sqrt( mean(depochs.^2,1) ) );
end
%{
This code performs a loop over a range of scales and computes the root-mean-square (RMS) over different time scales. 
Specifically, for each scale, it divides the data into non-overlapping epochs of that scale and computes the RMS over each epoch. 
The number of epochs is determined based on the length of the data and the scale. The detrend function is used to remove any linear 
trend within each epoch before computing the RMS. The resulting RMS values are stored in an array called rmses.
%}

% fit a linear model to quantify scaling exponent
A = [ ones(length(scales),1) log10(scales)' ];
dfa = (A'*A) \ (A'*log10(rmses)');
%{
This code performs a linear regression to fit a straight line to the log-transformed data, where the log-transformed RMS 
values are plotted against the log-transformed scales. The scaling exponent, which represents the slope of this line, can be 
interpreted as a measure of the long-range correlations in the data.

Taking the logarithm of the x-axis or y-axis (or both) is often done in power-law analysis to transform the data into a linear form. 
This transformation makes it easier to fit a linear model and estimate the slope of the line in the log-log space, which corresponds 
to the scaling exponent of the power-law relationship between the variables. This is useful because in the linear form the slope can be 
interpreted directly as the scaling exponent, whereas in the original non-linear form, the scaling exponent is not as straightforward to determine.
%}

% plot the 'linear' fit (in log-log space)
figure(6), clf
plot(log10(scales/srate),log10(rmses),'ks','linew',3,'markerfacecolor','w','markersize',18)

hold on
plot(log10(scales/srate),dfa(1)+dfa(2)*log10(scales),'k--','linew',2)

%% now repeat the procedure using DMA
%{
The log transformation is used to scale the data and make it more suitable for analysis. By taking the log of the data, 
we compress the range of values and reduce the influence of extreme values, making it easier to compare values 
across different scales. This is particularly useful for the analysis of time series data, where the range of values can be 
quite large and change over time. The log transformation is also useful for linearizing relationships between variables, 
making it easier to fit linear models and perform statistical tests.
%}

for scalei = 1:length(scales)
    
    % create kernel for this scale
    nConv  = N+scales(scalei)-1;
    kernel = fft(ones(scales(scalei),1)/scales(scalei), nConv);
    hfKrn  = floor(scales(scalei)/2)+1; % half of kernel
    %{
    nConv is the size of the padded convolution output, which is the sum of the signal length and the kernel length minus one. 
This ensures that the output is of the same size as the input signal.

kernel is a normalized boxcar filter of length scales(scalei) that is convolved with the signal. 
The boxcar filter is created by generating a vector of ones of length scales(scalei) and dividing it by the length of 
the filter scales(scalei). The division normalizes the filter, ensuring that the amplitude of the filter doesn't change as the length of the filter changes.

hfKrn is the index of the center of the kernel. It is used to ensure that the convolution output is of the same length 
as the input signal by discarding hfKrn samples from the beginning and hfKrn-1 samples from the end of the output signal.
    %}
    
    % mean-smooth as convolution
    convres = ifft( fft(perfdat, nConv) .* kernel );
    convres = convres(hfKrn:end-hfKrn+1+mod(nConv,2));
    %{
    In this code block, convres is computed using the convolution theorem. First, the input signal perfdat is Fourier transformed 
using fft, and a kernel of length scales(scalei) is created as a series of delta functions with amplitude 1 divided by the length of the kernel. 
The kernel is also Fourier transformed using fft. The product of the Fourier transformed perfdat and the Fourier transformed kernel is then 
taken element-wise using the .* operator. Finally, convres is obtained by taking the inverse Fourier transform of the product using ifft, and then 
removing the first and last hfKrn-1 elements of the result, where hfKrn is half the length of the kernel, to ensure that convres has the same length 
as the original signal perfdat.
    %}
    residX = perfdat - convres;
        
    n = floor(N/scales(scalei)); % number of epochs
    epochs = reshape(residX(1:n*scales(scalei)),scales(scalei),n)';
    rmses(scalei) = mean( sqrt( mean(epochs.^2,1) ) );
end

dma = (A'*A) \ (A'*log10(rmses)');


h1 = plot(log10(scales/srate),log10(rmses),'kd','linew',3,'markerfacecolor','w','markersize',18);

hold on
h2 = plot(log10(scales/srate),dma(1)+dma(2)*log10(scales),'k--','linew',2);

legend({'Data (DFA)';[ 'Fit (DFA m=' num2str(dfa(2)) ')' ];...
        'Data (DMA)';[ 'Fit (DMA m=' num2str(dma(2)) ')' ]})
xlabel('log10(Scales)')
ylabel('log10(RMS)')

set(h1,'color',ones(1,3)*.6)
set(h2,'color',ones(1,3)*.6)
set(gca,'xlim',[-.1 2.8],'ylim',[1.2 4.1])

%% local vs. global maximum using max

% create signal
perfdat = -40:.01:120;
signal = 3*sin(perfdat)./perfdat + sin(perfdat-20)./(perfdat-20);

figure(7), clf
plot(perfdat,signal,'k.-','linew',2)
xlabel('X'), ylabel('Y')
title('Two sincs, in synch')

% global maximum
[maxGval,maxGidx] = max(signal);
hold on
plot(perfdat(maxGidx),maxGval,'kv','linew',3,'markersize',15,'markerfacecolor','r')

% local maximum in specified range
range4max = [10 40];

rangeidx = dsearchn(perfdat',range4max');

[maxLval,maxLidx] = max(signal(rangeidx(1):rangeidx(2)));
% is this accurate??
plot(perfdat(maxLidx+rangeidx(1)),maxLval,'kv','linew',3,'markersize',15,'markerfacecolor','g')
zoom on

%% local minima/maxima

% find local maxima
peeks = find(diff(sign(diff(signal)))<0)+1;

[~,idx] = sort(signal(peeks));
peeks(idx(end-1:end)) = [];

plot(perfdat(peeks),signal(peeks),'ro','linew',2,'markersize',10,'markerfacecolor','y')

%% end.

