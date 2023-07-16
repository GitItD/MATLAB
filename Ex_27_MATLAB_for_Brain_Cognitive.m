%% 27.5.13
%{
The Fourier coefficients contain information about power and phase.
WeÅfve worked only with power in this chapter, but the phases are
important for images, as you will see in the next few exercises. Write
code to separate the phases from the magnitudes of all Fourier coefficients,
then reassemble them in the same order. Use the inverse Fourier
transform to reconstruct the Saturn image (use only one color dimension
to have a 2D image). ItÅfs not as simple as separating the real and
imaginary parts, because phase and magnitude both use real and imaginary
information. (Hint: A Fourier coefficient can be represented as
aeip, where a is the magnitude and p is the phase.) Make sure your
reconstruction is identical to the original image.
%}
% Load data
saturn = double( imread('saturn.png') )./255;
figure(1), clf
subplot(1,2,1); imagesc(saturn)
% Seperate phases from magnitudes of all Fourier... The magnitude Mf
% and phase Pf can be obtained from the real and imaginary parts of Xf
% using the following equations:
Xf = fft2(squeeze(saturn(:,:,1)));   % Xf = Mf * exp(i*Pf) is the fourier transform of signal xt
Mf =sqrt( real(Xf).^2 + imag(Xf).^2 );     % Magnitude part
Pf = atan2( imag(Xf), real(Xf) );       % Phase part
% Reassemble them in the same order through iFFT (use only one color
% dimension for 2D image). 
Pfi = ifft2(Pf);
Mfi = ifft2(Mf);
saturn_new = abs(ifft2(Xf.*(Pfi + Mfi))); 
% 
subplot(1,2,2); imagesc(saturn_new)
% End
%% 27.5.13 ChatGPT 
% Load Saturn image
saturn = double(imread('saturn.png'))/255;
% Compute 2D FFT of the image
saturnF = fft2(saturn);
% Get magnitudes and phases of all Fourier coefficients
Mf = abs(saturnF);              % Magnitude part
Pf = angle(saturnF);            % Phase part
% Reassemble magnitudes and phases in the same order
saturnF_recon = Mf .* exp(1i * Pf);
% Compute the inverse Fourier transform to reconstruct the image
saturn_recon = real(ifft2(saturnF_recon));
% Display the original and reconstructed images side by side
figure(2), clf
subplot(1,2,1);
imshow(saturn);
title('Original Image');
subplot(1,2,2);
imshow(saturn_recon);
title('Reconstructed Image');
%% 27.5.14
%{
Now that you know how to extract and recombine the phases, itÅfs time
to scramble them. Replace the phase angles from all coefficients with
random numbers. What is the appropriate range of random numbers
to use? Is the reconstructed image recognizable as Saturn?
%}
rnum = randn(size(saturnF)) * 0.1;
Pf = angle(saturnF) .* rnum;            % Phase part
% Reassemble magnitudes and phases in the same order
saturnF_recon = Mf .* exp(1i * Pf);
% Compute the inverse Fourier transform to reconstruct the image
saturn_recon = real(ifft2(saturnF_recon));
% Display the original and reconstructed images side by side
figure(3), clf
subplot(1,2,1);
imshow(saturn);
title('Original Image');
subplot(1,2,2);
imshow(saturn_recon);
title('Reconstructed Image');
% End
%% 27.5.14 ChatGPT
% Load data
saturn = double( imread('saturn.png') )./255;
% Get Fourier coefficients
Xf = fft2(saturn);
% Separate Magnitude and Phase
Mf = abs(Xf);
Pf = angle(Xf);
% Scramble phases
Pf_scrambled = rand(size(Pf)) * 2*pi - pi; % random values between -pi and pi
% Reassemble coefficients with scrambled phases
Xf_scrambled = Mf .* exp(1i * Pf_scrambled);
% Reconstruct image with scrambled phases
saturn_scrambled = real(ifft2(Xf_scrambled));
% Display original and scrambled images side by side
figure(4), clf
subplot(1, 2, 1)
imshow(saturn)
title('Original')
subplot(1, 2, 2)
imshow(saturn_scrambled)
title('Scrambled')
%% 27.5.15
%{
Now keep all of the original phases, but shift them all by some constant
offset, say, ÉŒ/3. You will need to write some code to make sure that the
angles remain in an expected range (hint: mod).
%}
% Create constant
modulus = mod(pi/3, 0);
% Shift
Pf = angle(saturnF) + modulus;
% Reassemble magnitudes and phases in the same order
saturnF_recon = Mf .* exp(1i * Pf);
% Compute the inverse Fourier transform to reconstruct the image
saturn_recon = real(ifft2(saturnF_recon));
% Display the original and reconstructed images side by side
figure(5), clf
subplot(1,2,1);
imshow(saturn);
title('Original Image');
subplot(1,2,2);
imshow(saturn_recon);
title('Reconstructed Image');
% End
%% 27.5.15 ChatGPT
% Seperate phases from magnitudes of all Fourier coefficients
Xf = fft2(saturn); 
Mf = abs(Xf);
Pf = angle(Xf);
% Shift all phases by pi/3
Pf_shifted = mod(Pf + pi/3, 2*pi);
% Reassemble the Fourier coefficients with the original magnitudes and shifted phases
Xf_new = Mf .* exp(1i * Pf_shifted);
saturn_new = ifft2(Xf_new);
% Display the new image
figure(6), imshow(saturn_new);
%% 27.5.16
%{
Finally, randomize only some phases while leaving others intact. Using
the box-shaped filter from the end of this chapter, randomize only the
phases inside that box. Is the reconstructed image recognizable? How
about when you randomize only the phases outside that box?
%}
% Get Fourier coefficients
Xf = fft2(squeeze( saturn(:,:,1) ));
% Separate Magnitude and Phase
Mf = abs(Xf);
Pf = angle(Xf);
% Scramble phases
Pf_scrambled = rand(size(Pf)) * 2*pi - pi; % random values between -pi and pi
% Reassemble coefficients with scrambled phases
Xf_scrambled = Mf .* exp(1i * Pf_scrambled);
% Reconstruct image with scrambled phases
saturn_scrambled = real(ifft2(Xf_scrambled));

% get sizes of image and midpoints
imgdims = size(Xf);
midX = round(imgdims(2)/2);
midY = round(imgdims(1)/2);
% size of the filter in pixels.
nPix2use = 100;
% create low-pass filter kernel
loPass2d = zeros(imgdims(1:2));
loPass2d(midY-nPix2use:midY+nPix2use, midX-nPix2use:midX+nPix2use) = 1;
% create high-pass filter kernel
hiPass2d = ones(imgdims(1:2));
hiPass2d(midY-nPix2use:midY+nPix2use, midX-nPix2use:midX+nPix2use) = 0;
% original image
figure(6), clf, subplot(131)
imagesc(saturn)
axis off, axis square
title('Original')
% inside the box
filtimg = real(ifft2( Xf_scrambled.*fftshift(loPass2d) ));
subplot(132)
imagesc(filtimg)
axis off, axis square
title('Scrambled In-box')
% outside the box
filtimg = real(ifft2( Xf.*fftshift(hiPass2d) ));    % This needs to be corrected, since you're not using high-pass filter, but outside the box filter
subplot(133)
imagesc(filtimg)
axis off, axis square
set(gca,'clim',[-.2 .2]/10)
title('Scrambled Out-box')