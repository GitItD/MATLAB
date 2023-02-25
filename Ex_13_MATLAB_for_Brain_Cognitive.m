clear
clc
%{
The codes below were created by an AI (ChatbotGPT) to test its capabilities
in programming. While it was able to produce certain codes that did work,
it also struggled to answer the questions correctly, and may other codes did not
work. For general educational purposes, ChatbotGPT may be helpful for
learning to program, but for specific questions, close monitoring should be
used to make sure the AI performs up to standard. 

More future analysis may be needed to determine whether if ChatbotGPT can
be used for easing of coding in MATLAB. 

One way I have found ChatbotGPT to be useful is in helping to simpify
learning the codes. At certain times, it is difficult to analyze codes and
to get an intuitive feelings for what they mean. ChatbotGPT had been able
to translate the codes into simplified English, which allowed me to
accelerate my learning.
%} 

%% 13.6.1:
%{ 
Generate a time series of random numbers, and up-sample the time
series using two methods of interpolation: the zero-padding theorem
and griddedInterpolant. Make sure the number of time points is
the same. How comparable are the results, and is it possible to find
parameters to make the results more similar?

To generate a time series of random numbers, we can use the built-in Matlab function "rand" 
which generates a random number between 0 and 1. To create a time series, 
we can multiply the random number by a value, such as 100, 
and add a value such as 50 to create a time series with values between 50 and 150.
For example:
%}
times = 1:10; original_data = 50 + (100*rand(1,10));
%{
To up-sample the time series using zero-padding, we can use the "interp" function in Matlab, 
which allows us to specify the number of points we want to interpolate the data to. 
In this case, we want to increase the number of time points by a factor of 10, so 
we can set the number of points to be 10 times the original number of points.
For example:
%}
new_times = 1:0.1:10; zero_padded_data = interp(original_data,10);
%{
To up-sample the time series using griddedInterpolant, we can use the "griddedInterpolant" function in Matlab, 
which allows us to specify the interpolation method. In this case, we will use linear interpolation.
For example:
%}
F = griddedInterpolant(times,original_data,'linear'); gridded_data = F(new_times);
%{
The results of the two interpolation methods may not be exactly the same, as the zero-padding theorem 
and linear interpolation use different methods to estimate the values at the new time points. 
However, by adjusting the interpolation method in griddedInterpolant, it may be possible to make the results more similar.
 For example, using a more advanced interpolation method such as spline interpolation 
may produce results that are more similar to the zero-padding theorem.
Additionally, it is possible to fine-tune the parameters of the interpolation 
methods, such as the order of the polynomial used in interpolation, or the 
smoothing parameter used in spline interpolation, to make the results more similar.
However, this would require further experimentation and analysis of the data.
%}

%% 13.6.2
%{
Simulate 3 seconds of random data (use randn) at 1 kHz. Down-sample
the data to 200 Hz. Then up-sample the data back to 1 kHz using linear
and spline interpolation. Plot all of the time courses in the time domain
and their power spectra. What is the effect of down-sampling then upsampling
on the resulting time series and the power spectrum, and
does linear or spline interpolation more closely match the original time
series?
%}

% Simulating 3 seconds of random data at 1 kHz can be done using the following code:
fs = 1000; % Sampling frequency
t = 0:1/fs:3; % Time vector
originalData = randn(size(t)); % Generate random data

% To down-sample the data to 200 Hz, we can simply take every 5th sample from the original time series:
downSampledData = originalData(1:5:end);

% To up-sample the data back to 1 kHz, we can use the interp1 function to perform linear and spline interpolation:
% Linear interpolation
fsUp = 1000; % Upsampled frequency
tUp = 0:1/fsUp:3; % Upsampled time vector
linearInterpData = interp1(1:5:length(originalData), downSampledData, 1:length(tUp), 'linear');

% Spline interpolation
splineInterpData = interp1(1:5:length(originalData), downSampledData, 1:length(tUp), 'spline');

% To plot the time courses in the time domain, we can use the following code:
figure(1);
plot(t, originalData, 'k');
hold on;
plot(tUp, linearInterpData, 'b');
plot(tUp, splineInterpData, 'r');
legend('Original', 'Linear Interpolation', 'Spline Interpolation');
xlabel('Time (s)');
ylabel('Amplitude');

% To plot the power spectra, we can use the following code:
figure(2);
originalDataFft = abs(fft(originalData));
linearInterpDataFft = abs(fft(linearInterpData));
splineInterpDataFft = abs(fft(splineInterpData));

plot(linspace(0, fs, length(originalDataFft)), originalDataFft, 'k');
hold on;
plot(linspace(0, fs, length(linearInterpDataFft)), linearInterpDataFft, 'b');
plot(linspace(0, fs, length(splineInterpDataFft)), splineInterpDataFft, 'r');
legend('Original', 'Linear Interpolation', 'Spline Interpolation');
xlabel('Frequency (Hz)');
ylabel('Power');

%{
The effect of down-sampling then up-sampling on the resulting time series is that it 
introduces distortion and noise. The power spectrum is also affected, and the original 
high-frequency components are lost due to the down-sampling process. 
The linear interpolation method tends to introduce more distortion and noise compared 
to the spline interpolation method, resulting in a less accurate representation of the original
time series. The spline interpolation method tends to introduce less distortion and noise, 
resulting in a more accurate representation of the original time series. It's possible to find 
parameters to make the results more similar by using a higher order interpolation 
method or using a more sophisticated interpolation algorithm.
%}