%% Matlab for Brain and Cognitive Scientists
% Matlab code for Chapter 12
% Mike X Cohen
% 
% This code accompanies the book 
%      "Matlab for Brain and Cognitive Scientists" (MIT Press, 2017).
% Using the code without following the book may lead to confusion, 
% incorrect data analyses, and misinterpretations of results. 
% Mike X Cohen assumes no responsibility for incorrect use of this code. 

%{
Study the codes below after you understand at minimum, the fundamentals of 
convolution and fast fourier transform.
%}
clear
clc

%% "manual" time-domain convolution

% signal
signal = zeros(1,30);
signal(11:19) = 1;

% kernel (Gaussian)
kernel = [0, 0.31, 0.77, 1, 0.77, 0.31, 0];

% normalize the kernel
kernel = kernel./sum(kernel);


% sizes
N = length(signal);               % length time of the signal
M = length(kernel);             % length time of the kernel
halfKern = floor(M/2);       % Nyquist limit sampling rate, floor rounds M/2 toward negative infinity (ex: 1.9 = -2.000, 3.8 = 3)

%{
Data that we'll use for convolution (must be zero-padded).
In order to make the dot product a valid computation, the signal needs
to be zero-padded at the beginning by M ? 1 points, and at the end by M ? 1
points (minus 1 because the kernel and signal overlap by one point. 
The resulting dot product from this step of convolution
is plotted in a location corresponding to the center of the kernel
relative to where it is lined up with the signal. For this reason, it is useful to
make sure the kernel has an odd number of points.
%}
dat4conv = [zeros(1, M-1), signal];

% initialize convolution output
conv_res = zeros(1, N+M-1);

%{ 
run convolution
The rest of convolution is simple. At each step of convolution, keep the
signal fixed and slide the kernel one time step to the right. Repeat the dot
product computation. Slide, dot product, repeat. Over and over again, until
the left-most point of the kernel is aligned with the right-most point of the
data.
%}
for ti = M+1:N-M+1
    tempdata = dat4conv(ti:ti+M-1);
    
    % compute dot-product (don't forget to flip the kernel backwards!)
    conv_res(ti) = sum( tempdata.*kernel(end:-1:1) );
end

%{ 
trim the wings
The result of convolution starts before the signal begins and ends after
the signal ends, meaning that the result of convolution is longer than the
original signal. I call these the ÅgwingsÅh of convolution. How long is the
result of convolution? The result of convolution is half the length of
the kernel too long at the beginning and half the length of the kernel too
long at the end. You might initially think the result is length N + M, but itÅfs
actually N + M ? 1. To understand why the Ågminus 1Åh is necessary, itÅfs easiest
to see it visually explained, which I hope is clear enough in figure 12.2.
N + M ? 1 is an important quantity, and you should not forget it.
Technically, convolution is now finished. However, it is useful to cut off
the wings for simplicity.

If the convolution is not done with N + M ? 1 points, the convolution
wings wrap around such that the end of the time series sums into the beginning
of the time series, and vice versa. This is called circular convolution
(using N + M ? 1 is called linear convolution). If the time series is long relative
to the length of the kernel, this awkward wraparound effect has minimal
impact except at the edges. YouÅfll see an example of circular convolution
in chapter 20.
%}
conv_res = conv_res(halfKern+1:end-halfKern);

figure(1), clf

% plot signal
subplot(311)
plot(signal)
set(gca,'xlim',[0 N],'ylim',[min(signal)-.05 max(signal)+.05])

% plot kernel
subplot(312)
plot(kernel)
set(gca,'xlim',[0 N])

% plot convolution result
subplot(313)
plot(conv_res,'o-','linew',2,'markerface','g','markersize',9)
set(gca,'xlim',[0 N],'ylim',[min(signal)-.05 max(signal)+.05])

% Run this next line to show that our manual convolution 
% matches the output of Matlab's conv function.
%hold on, plot(conv(signal,kernel,'same'),'r*')

%% time-domain convolution as frequency-domain multiplication
%{
Time-domain convolution is okay, but (1) itÅfs slow, and (2) itÅfs not very
intuitive why it is so powerful as a signal processing technique (at least, not
to me). Convolution can also be implemented in the frequency domain,
and this solves both problems. Before learning how to implement convolution
in the frequency domain, you need to know about the convolution
theorem.

The convolution theorem states that convolution in the time domain
is equivalent to multiplication in the frequency domain. This means that
we can implement convolution in two distinct ways and arrive at the
same result (figure 12.3). The full proof of the convolution theorem is
not presented here, but briefly, it involves demonstrating that the
Fourier transform of a convolution between signals X and Y can be
reduced to the multiplication of the Fourier transform of X by the Fourier
transform of Y.

1. Time domain: Compute repeated dot products between the kernel and
temporally corresponding points in the signal. At each repeat, shift the
kernel over one time step relative to the signal.

2. Frequency domain: Compute the FFT of the signal, the FFT of the
kernel, point-wise (i.e., frequency by frequency) multiply the Fourier
coefficients together, and compute the inverse FFT.

Important! No matter which way you perform convolution, the result of
convolution must have length N + M ? 1. If you perform convolution in the
frequency domain, you need to set the length of the FFT to be N + M ? 1,
but you do not need to zero-pad the start and end of the signal the way you
would in the time domain.
%}

% Our data (same as before, incase you forgot)
signal = zeros(1,30);
signal(11:19) = 1;
kernel = [0, 0.31, 0.77, 1, 0.77, 0.31, 0];

% step 1: get the right sizes
nData = length(signal);             % length of signal N
nKern = length(kernel);            % length of kernel M
nConv = nData+nKern-1;         % length of convolution result = N + M - 1
nHfkn = floor(nKern/2);           % half the length of kernel, the Nyquist sampling rate

% step 2: FFTs
signalX = fft(signal,nConv);        % FFT of signal
kernelX = fft(kernel,nConv);      % FFT of kernel
hold on, plot(conv(signal,kernel,'same'),'r*')         

% step 3: normalize the kernel
%{
Depending on the kernel, it may require normalization in the frequency
domain to ensure that the result of convolution will be in the
same units as the original signal. For wavelets (see chapter 19), this can
be implemented by dividing the Fourier spectrum by its maximum.
%}
kernelX = kernelX./max(kernelX);

% step 4: point-wise multiplication
%{
Point-wise multiply the two Fourier spectra and take the inverse FFT.
This point-wise multiplication will fail unless both FFTs are the same
length, which is a good check that you did step 2 correctly.
%}
convres = ifft( signalX.*kernelX );

% step 5: trim wings
%{
Point-wise multiply the two Fourier spectra and take the inverse FFT.
This point-wise multiplication will fail unless both FFTs are the same
length, which is a good check that you did step 2 correctly.
%}
convres = convres(nHfkn+1:end-nHfkn);

hold on
plot(convres,'rd','linew',4)

%{
Discussion on the result:

Steps 3 and 5 are not technically necessary. Step 3 is useful if you want
to keep the postconvolution results in the same scale as the original data
(e.g., É V for EEG data). Concerning step 5, it is difficult to imagine a
situation where you would want to keep the convolution wings, and leaving
them attached can only increase the possibility of confusion and
mistakes.

In step 4, thinking about convolution as point-wise multiplication between the
frequency spectra of the signal and that of the kernel offers an enlightening
perspective on what convolution does and why it works. LetÅfs start
with a simple example of random noise (the signal) and a Gaussian (the
kernel).

The figure plots show the two time series and the result of their convolution.
The Gaussian acts as a low-pass filter. This happens because the shape of the
Gaussian in the frequency domain is a negative exponential, and when the
two spectra are point-wise multiplied, only the frequencies in the signal
that match those with non-zero power in the Gaussian are preserved in the
inverse FFT. More specifically, the higher frequencies are attenuated because
the power spectrum of the Gaussian is nearly zero. YouÅfll see additional
examples of convolution in chapter 20. The top row shows the time domain,
and the bottom row shows the frequency domain. The left-most panels show the
time series signal; the middle panels show the kernel; and the right-most panels
show the result of convolution. Notice that the result of convolution is a mixture of
the characteristics of the signal and the kernel.
%}

%% convolution between signal and noise

% signal of random noise
signal = randn(1,1000);

% gaussian kernel (width is hard-coded)
x = -2:.01:2;
kernel = exp( (-x.^2) / .1 );


% step 1: get the right sizes
nData = length(signal);
nKern = length(kernel);
nConv = nData+nKern-1;
nHfkn = floor(nKern/2);

% step 2: FFTs
signalX = fft(signal,nConv);
kernelX = fft(kernel,nConv);

% step 3: normalize the kernel
kernelX = kernelX./max(kernelX);

% step 4: point-wise multiplication
convres = ifft( signalX.*kernelX );

% step 5: trim wings
convres = convres(nHfkn+1:end-nHfkn);

figure(6), clf
subplot(231)
plot(signal)
set(gca,'xlim',[0 nData])
xlabel('Time (a.u.)')
title('Time domain: signal')

subplot(232)
plot(kernel)
set(gca,'xlim',[0 nData])
xlabel('Time (a.u.)')
title('Time domain: kernel')

subplot(233)
plot(convres)
set(gca,'xlim',[0 nData])
xlabel('Time (a.u.)')
title('Time domain: convolution result')


subplot(234)
plot(2*abs(signalX))
set(gca,'xlim',[0 nConv])
xlabel('Frequency (a.u.)')
title('Frequency domain: signal')

subplot(235)
plot(2*abs(kernelX))
set(gca,'xlim',[0 80])
xlabel('Frequency (a.u.)')
title('Frequency domain: kernel')

subplot(236)
plot(2*abs(signalX.*kernelX))
set(gca,'xlim',[0 80])
xlabel('Frequency (a.u.)')
title('Frequency domain: convolution result')

%% 2D convolution
%{
Convolution in two dimensions is used in image processing. The principle
is the same as convolution in one dimension, except that it is done over
two dimensions. As with 1D convolution, 2D convolution can be done the
slow way (e.g., the space domain) using two nested for-loops, or it can be
done the fast and efficient way using frequency-domain multiplication
(make sure to use fft2, not fft!). IÅfll just focus on the frequency-domain
implementation here.

Our low-pass filter will be a 2D Gaussian. Although we have an extra dimension 
to deal with, the basic five-step plan of convolution is the same. 
There are a few minor departures from the 1D process:

From Step 1:                nConv = N+M-1;
From Step 2:                picX = fft2(pic,nConv(1),nConv(2));
From Step 5:                cr=cr(halfK(1)+1:end-halfK(1),halfK(2)+1:end-halfK(2));

N and M are now two-element vectors (height and width of the image)
instead of one-element scalars. In the fft2 function, there are two inputs
to specify the nfft parameter?one for height and one for width.
And be careful that the wings must be trimmed on both dimensions. The
variable cr is a 2D matrix, so there will be 2 indexes instead of 1.
%}

% import picture (We'll use only the red layer, because today feels like a red day, don't you think?)
pic = imread('saturn.png');
pic = double(squeeze(pic(:,:,1)));

% create gaussian
[x,y]  = meshgrid(-250:250);
s = 30;
gaus2d = exp( -(x.^2 + y.^2)/(2*s^2) );


% step 1
N     = size(pic);
M     = size(gaus2d);
nConv = N+M-1;
halfK = floor(M/2);

% step 2
picX  = fft2(pic,nConv(1),nConv(2));
gausX = fft2(gaus2d,nConv(1),nConv(2));

% step 3
gausX = gausX ./ max(gausX(:));

% step 4
cr = ifft2( picX.*gausX );

% step 5
cr = cr( halfK(1)+1:end-halfK(1) , halfK(2)+1:end-halfK(2) );


figure(2), clf
subplot(231)
imagesc(pic)
axis square, axis off
title('Saturn in space')

subplot(232)
imagesc(gaus2d)
axis square, axis off
title('Gaussian in space')

subplot(233)
imagesc(cr)
axis square, axis off
title('Convolution result in space')


subplot(234)
imagesc(fftshift(log10(abs(picX))))
set(gca,'clim',[2.5 5])
axis square, axis off
title('Saturn in frequency')

subplot(235)
imagesc(fftshift(log10(abs(gausX))))
set(gca,'clim',[-20 0])
axis square, axis off
title('Gaussian in frequency')

subplot(236)
imagesc(fftshift(log10(abs( picX.*gausX ))))
set(gca,'clim',[-20 0])
axis square, axis off
title('Convolution result in frequency')

%% end.
