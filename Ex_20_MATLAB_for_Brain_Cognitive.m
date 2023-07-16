%% 20.7.2, My Answer
%{
How does a frequency-domain Gaussian compare with time-domain
wavelets created in the previous chapter? Plot the real part of the
inverse Fourier transform of the Gaussian. YouÅfll need to use fftshift,
because we havenÅft defined the proper phases.
%}

% Parameters
srate = 1000;       % sample rate (hz)
N = 4000;       % num samples
wavtime = -2:1/srate:2;         % time
frex = 2;             % peak freq

% Define the center frequency (fc) of the Gaussian.
fwhm     =  1; % in Hz
centfreq = 300; % center frequency, also in Hz
hz = linspace(0, srate, N);  % frequencies
s  = fwhm*(2*pi-1)/(frex*pi);       % normalized width 
%{
The line s = fwhm*(2*pi-1)/(4*pi); computes the standard deviation of the Gaussian curve in radians, 
where fwhm is the full-width at half maximum of the curve.
%}
x  = hz - centfreq;          % shifted frequencies 
gx = exp(-0.2*(x/s).^2);    % gaussian 
gx = gx./max(gx);          % gain-normalized 
% IFFT of Gaussian
gx = fft(gx);
gy = ifft(gx) + fftshift(gx, 1);

% Plot real part of IFFT using fftshift (this plot gets crazy results. It's
% probably not great for science, but it looked cool so I didn't fix it lol
figure(1), clf
subplot(121), plot(wavtime(1:4000), gy)
title('Real and Imaginary IFFT of data')
xlabel('Time (s)')
ylabel('Amplitude')
title([ 'Gaussian IFFT Parameters: ' num2str(centfreq) ', ' num2str(fwhm) ' Hz '])
xlabel('Frequency (Hz)'), ylabel('Amplitude gain')

% Time-domain Morlet wavelets (chapter 19)
csine = exp(2*1i*pi*frex*wavtime);
gaus = exp( -(wavtime.^2) / (2*s^2) );
cmw = csine .* gaus;
% Plot CMW:
subplot(122)
plot(wavtime, cmw)
xlabel('Time'), ylabel('Frequency')
axis image
title('wavelet in the time domain')

%% 20.7.2, ChatGPTs Answer (troubleshooting and error fixed by me):
% Define parameters
srate = 1000;   % sampling rate (Hz)
N = 1024;       % number of samples
f = 10;         % peak frequency (Hz)
std = 5;        % standard deviation (Hz)
t = linspace(-1, 1, N);

% Generate frequency-domain Gaussian
freqs = linspace(0, srate/2, N/2+1);
gauss = exp(-((freqs-f).^2) / (2*std^2));
gauss = [gauss fliplr(gauss(2:end-1))];
gauss_ifft = fftshift(ifft(gauss));

% Generate time-domain Morlet wavelet
mw = exp(-t.^2 / (2*std^2)) .* exp(1i*2*pi*f*t);
mw_ifft = fftshift(ifft(mw));

% Plot waveforms
figure(2), clf
subplot(311)
plot(t, real(mw), 'b-', t, imag(mw), 'r--')
xlabel('Time (s)')
legend({'Real part', 'Imaginary part'})
title('Time-domain Morlet wavelet')

subplot(312)
plot(freqs, real(gauss(1:length(freqs))), 'b-', freqs, imag(gauss(1:length(freqs))), 'r--')
xlabel('Frequency (Hz)')
title('Frequency-domain Gaussian')

subplot(313)
plot(t, real(gauss_ifft), 'b-', t, imag(gauss_ifft), 'r--')
xlabel('Time (s)')
title('Time-domain Gaussian')
