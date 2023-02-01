%% Fourier transform examples
%{
1) Construct a complex sine wave A
2) Compute the complex dot product between that sine wave A and the signal (compute the Fourier Series Coefficient)
3) Construct another complex sine wave B at a different frequency 
4) Compute the complex dot product of complex sine wave B
5) This process repeats for as many frequencies as there are time points in the data. The number of frequencies in the discrete Fourier transform is
defined by the number of time points.
%}

L = length(signal);     % The length of signal, our source data
fTime = (0:L - 1)/L;    % The length of time

for fr = 1:L    % The length of the frequency is the length of the signal
    fSine = exp(-i*2*pi*(fr - 1).*fTime); % Computation of the Fourier Series Coefficients
    signalX(fr) = sum(fSine.*signal);   % Integral summation over time to complete the series
end
signalX = signalX/L;        % Dividing the total sum with scaling factor matches the Fourier coefficient to the original signal. 

% Nyquist limit: 
%{
After computing the Fourier transform, you will want to label the Fourier
coefficients in terms of frequencies in Hz.  
%}
nyquist = signalX/2;      % Base of Nyquist limit is signal rate divided by 2.
hz = linspace(0, nyquist, L/2+1);       % So the Hertz limit is from 0 to Nyquist, and should be intervaled at rate greater than N/2.

%% Fast Fourier transform
%{
It is important to understand how the discrete Fourier transform is implemented,
because thatÅfs the best way to understand how the Fourier transform
works. In practice, however, you should use the fast Fourier transform
(FFT), which can save your analysis code hours or even days of computation
time.

There are several algorithms for computing the fast Fourier transform, but
in general they work by putting all of the sine waves into one matrix rather
than using a for-loop, using matrix factorizations to break up that large
matrix into many simpler matrices with many zeros (sparse matrices), and
then reconstructing the Fourier spectrum.

In the example below, we'll use FFT to find the frequency component of a
signal buried in noise. The sampling frequency is 1 kHz and the signal
duration is 1.5 sec.
%} 

fr = 1000;      % sampling frequency of 1 kHz
sr = 1/fr;      % sampling period
L = 1500;       % if the signal duration is 1.5 s, then by definition of hz, 1 khz = 1.5 * 1000 = 1500
t = (0:L-1)*sr;      % time vector

a = 0.7;         % The amplitude of our complex sinewave equation
f = 50;          % The frequency in hz of our complex sinewave equation
f2 = 120;        % The 2nd frequency in hz of our complex sinewave equation
csw = a*sin(2*pi*f*t) + sin(2*pi*f2*t);     % Our complex sinewave series
s_X = csw + 2*randn(size(t));     % Add a noise signal to complex sinewave equation

figure
plot(fr*t(1:f), s_X(1:f))       % Plot the noisy signal in the time domain. 
title('Signal corrupted with zero-mean random noise')
xlabel('t (milliseconds)')
ylabel('X(t)')

Y = fft(s_X);       % compute the Fourier transform of the noisy signal
P2 = abs(Y/L);  % the absolute value of Fourier transform divided by signal duration is spectrum 2
P1= P2(1:L/2+1);        % computing spectrum 1 based on P2 and even-valued signal length L
P1(2:end - 1) = 2*P1(2:end - 1); 
%{ 
Define the frequency domain f and plot the single-isded amplitude spctrum
P1. The amplitudes are not exactly at 0.7 and 1 as expected, becaue of the
added noise. On average, longer signals produce better frequency
approximations.
%}
fs = fr*(0:(L/2))/L;

figure
plot(fs, P1)
title('Single-sided amplitude spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')

% Now take the Fourier transfom of the original, uncorrupted signal and retrieve the exact amplitudes, 0.7 and 1.0
Y = fft(csw);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

figure
plot(fs, P1)
title('Single-sided amplitude spectrum of csw(t)')
xlabel('fs (Hz)')
ylabel('|P1(fs)|')

%% 2D Fourier Transform
%{
It is possible to compute Fourier transforms of more than 1D by the
command, Y = fftN(x), where N is the number of dimensions of the data. In
the example below, we show how to compute a 2D Fourier transform. 
%}
pic = imread('saturn.png')

tic
picX = fft2(pic);
a = toc
tic
picX1 = fft(pic);
picX1 = fft(picX1, [], 2);
b = toc

tf = sum(picX(:)) == sum(picX1(:))

%{
The variables picX and picX1 are the same, although MATLABÅfs fftn
function is a bit faster than calling the fft twice. The resulting Fourier coefficients
come in a 2D matrix. ItÅfs a bit confusing at first to interpret the
results, but keep in mind that the power spectrum shows the power along
each dimension separately. For this reason, it is useful to think about the 2D
power spectrum not as going from left to right, but as going from the corners
of the matrix to the center of the matrix. And for this reason, it is common
to shift the coefficients around so the four corners (corresponding
to low spatial frequencies) meet at the middle, and the original center 
(corresponding to high spatial frequencies) is shifted to the corners. This
shifting is so common that MATLAB has a built-in function to accomplish
this, called fftshift.
%}