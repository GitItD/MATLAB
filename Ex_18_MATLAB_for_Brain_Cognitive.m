%% 18.9.5
%{
Figure 18.1 showed the effect of an edge in the time domain on the
power spectrum. Is this detrimental for the spectrum of sinusoidal
components? To find out, add a pure sine wave to the variable ts.
%}
% create time series with sine waves
ts = zeros(100,1);
time = 0:1/(52-48):1
ts(48:52) = 0.5*sin(2*pi*10*time + 10);
% create time series with edges (for comparison)
ts2 = zeros(100,1);
ts2(48:52) = 10;

figure(1), clf

% plot time-domain signal
subplot(211)
hold on
plot(ts,'k','linew',3)
xlabel('Time (a.u.)'), ylabel('Amplitude (a.u.)')
set(gca,'ylim',[-.5 10.5])
plot(ts2,'--','linew',3, 'markerface', 'g')
legend('Pure Sinewave', 'Time Series + Edges')
hold off

% plot power spectrum of signal
subplot(212)
hold on
plot(abs(fft(ts)/100).^2,'ro','linew',1.5,'markerface','g')
plot(abs(fft(ts2)/100).^2, '--', 'LineWidth', 3);
xlabel('Frequencies (a.u.)'), ylabel('Power (a.u.)')
legend('Pure Sinewave', 'Time Series + Edges')
hold off
%% 18.9.6
%{
The online MATLAB material includes a file called mouseHippocampus.
mat, which contains 100 trials of LFP recordings from the hippocampus.
On each trial, a visual stimulus was displayed. Show the power
spectrum from this electrode for two methods of averaging over trials.
First, compute the FFT and extract the power spectrum of each trial,
and then average the power spectra together. Second, average the timedomain
LFP traces together, and then compute the power spectrum.
Show the results on the same plot.
%}
%% 18.9.7
%{
How much of the data are lost when tapering? To find out, create a 100-
by-120 matrix of normally distributed random numbers (120 epochs,
each with 100 time points). Then, taper each epoch with a Hann window
(can you do this without a loop?). Make sure you are tapering the
time dimension, not the epoch dimension. At each time point over
epochs, compute the sum of squared errors between the tapered signal
and the original signal. The result will be a 1-by-100 vector of how
much and where information is lost due to tapering. Finally, repeat this
procedure using Hamming and Gaussian windows. Plot the results for
different tapers overlaid on the same window.
%}
%% 18.9.9
%{
Compute the power spectrum of the online data EEGrestingState.mat.
The variable eegdata is a time-by-epoch matrix. Use the entire time
series in one FFT, with and without applying a Hann taper. In one figure
with 2 Å~ 1 subplots, plot in the upper subplot the time-domain
data before (black line) and after (red line on top) applying the Hann
taper. In the lower subplot, plot the power spectrum from the entire
time series without (black line) and with (red line on top) the Hann
taper.
%}
%% 18.9.10
%{
With the same resting-state data, cut the time series into nonoverlapping
epochs of two seconds. Then taper each epoch using a
Hann window, take the power spectrum, and then average the power
spectrum over epochs. First do this using a loop over epochs. Then get
rid of the loop by inputting a matrix into the fft function. Check that
the loop and no-loop code produce identical results.
%}
%% 18.9.11
%{
Using the epoched resting-state data, compute the power spectrum
using an N parameter for the FFT corresponding to the length of the
entire time series. That is, zero-pad the FFT of the epochs so the frequency
resolution is the same for the epochs and for the continuous
time series. Then plot three power spectra on the same plot (donÅft taper
the data for this exercise): from the continuous data, from the zeropadded
epoched data, and from the non-zero-padded epoched data.
You might want to scale up the power from the continuous data in
order to make the power spectra more directly comparable.
%}
%%18.9.15
%{
One way to attenuate the 1/f characteristic of the power spectrum is to
take the FFT of the first derivative of the time series (sometimes called
Ågpre-whiteningÅh). The MATLAB function diff takes the derivative,
and the result is length N ? 1 (this makes sense: the discrete derivative
of [2 1 4] is [?1 3]). In this problem, you will simulate ÅgpinkÅh noise with
1/f characteristics, add some sine waves, and then compute the power
spectrum before and after pre-whitening. One method to create 1/f
noise is to modulate the power spectrum of white noise and then take
its inverse Fourier transform. The code below does most of the work
for you.
%}
%%18.9.18
%{
What is the error in the following code? Write code to perform the
following sanity check: Create a 1-second sine wave at 13 Hz using a
sampling rate of 1 kHz. Then plot the power spectrum of that sine wave
and identify the frequency of the peak power. If the peak is not 13 Hz,
then the frequencies vector is incorrect. You should use the datacursormode
tool to check the exact frequency

hz = linspace(0,srate/2,floor(n/2));
%}