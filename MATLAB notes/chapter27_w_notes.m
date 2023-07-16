%% Matlab for Brain and Cognitive Scientists
% Matlab code for Chapter 27
% Mike X Cohen
% 
% This code accompanies the book 
%      "Matlab for Brain and Cognitive Scientists" (MIT Press, 2017).
% Using the code without following the book may lead to confusion, 
% incorrect data analyses, and misinterpretations of results. 
% Mike X Cohen assumes no responsibility for incorrect use of this code. 

%{
Background info:
2D mean filtering is a type of image filtering that is used to reduce noise and smooth 
an image by averaging the pixel values in a local neighborhood around each pixel. 
In this process, each pixel is replaced with the mean value of its surrounding pixels in a 
defined neighborhood, typically a square or circular window of a specified size.

Mean filtering is often used to remove random noise from an image, such as Gaussian 
or salt-and-pepper noise, while preserving the sharpness and edges of the image. 
The amount of smoothing or blurring effect depends on the size of the neighborhood 
window used for the filter.

However, one drawback of mean filtering is that it can also blur out fine details 
and edges in an image, which may be undesirable in certain applications. In these 
cases, more advanced filtering techniques such as median filtering or Gaussian 
filtering may be more appropriate.

2D median filtering is a method of noise reduction in digital images. It involves 
replacing the value of each pixel with the median value of its neighboring pixels 
within a specified window size. Unlike mean filtering, which replaces the value of 
each pixel with the average of its neighbors, median filtering is less sensitive to extreme 
values (outliers) and preserves the edges and details of an image.

Median filtering is often used to remove salt-and-pepper noise, a type of image noise that 
appears as white and black dots or speckles in an image. It can also be used to reduce 
Gaussian noise, another common type of image noise caused by random variations in pixel values.

Gaussian kernel smoothing is a type of image filtering that is commonly used for smoothing 
images in order to remove noise and to enhance edges. The process involves convolving the 
image with a Gaussian kernel, which is a type of mathematical function that can be used to 
represent a normal distribution.

An isotropic Gaussian filter is a filter that has the same standard deviation in all directions. 
This type of filter is used when the noise in the image is isotropic, meaning that it has the 
same amount of noise in all directions.

An anisotropic Gaussian filter, on the other hand, is a filter that has different standard 
deviations in different directions. This type of filter is used when the noise in the image is 
anisotropic, meaning that it has different amounts of noise in different directions.
%}

%% load in the Saturn picture

saturn = double( imread('saturn.png') )./255;
% imagesc(saturn)

% first, add some noise
noiselevel = .5;
saturn = saturn + randn(size(saturn))*noiselevel;

% let's see how noisy it is
figure(1), clf
imagesc(saturn), axis square

%% mean-smoothing the wrong way
%{
A 2D mean-smoothing filter is the same concept as?and has a very 
similar implementation to?the 1D running-average filter you learned 
about in chapter 20. ItÅfs a form of low-pass filtering and thus produces 
an image that is smoother (blurrier) than the original image. Mean filtering is an oft-used
denoising method because many sources of noise produce high spatial frequency artifacts.

Each pixel in the image is replaced by average of surrounding K-nearest neighbors algorithm pixels. 
%}
% What are the problems with this code?

k=9;
for i=k+1:size(saturn,1)
    for j=k+1:size(saturn,2)
        temp = saturn(i-k:i+k,j-k:j+k);
        saturn(i,j) = mean(temp(:));
    end
end
%{
My guess to the problems:
1. This code is running i and j of Saturn's 1D and 2D data, then extracting
those data and generating a mean out of the extracted data. However,
because this is a loop within a loop, we are increasing the numbers of
loops by i*j, which do not allow an accurate calculation of the mean. 
2. This code generates mean at temp; which goes from i-k to i+k, and j-k to
j+k. This means that we are generating the mean of the image in both rows
and columns. However, we are averaging the entire image when we use temp(:)
given that both rows and columns are added together. This doesn't allow us
to average based on surrounding means. 
3. The size of saturn variable is assigned to the generation of random
noise. This means that the mean can fluctuate greatly by the noise, which
makes for a poor filter. 

The answer:
1) The loop should stop k units before the end of the image, or size(saturn, 1)-k, 
to avoid index exceeds matrix dimensions error.
2) The image of Saturn is a 3D matrix, but the code is indexing it like a 2D matrix, 
leading to undesired results. The solution is to add another loop over elements 
in the third dimension or to smooth only one dimension.
3) As the center of the filter moves around the image, parts of the image get 
double-filtered, leading to erroneous behavior. The solution is to filter the image 
into a different variable, like saturnfilt(i,j,1) = mean(temp(:)).
4) The code overwrites the image with a noisier version of itself, losing the original 
data. A better implementation would be to create a new variable and preserve the original.
%}
%% mean smoothing the right way

k=9;

saturn = double( imread('saturn.png') )./255;
saturnFilt = saturn;

for dimi=1%:size(saturn, 3)
    for i=k+1:size(saturn, 1)-k
        for j=k+1:size(saturn, 2)-k
            temp = saturn(i-k:i+k, j-k:j+k, dimi);
            saturnFilt(i, j, dimi) = mean(temp(:));
        end % end loop over rows
    end % end loop over columns
end % end loop of dimensions

figure(2), clf
subplot(121)
imagesc(squeeze(saturn(:,:,1)))
axis xy, axis square, axis off
set(gca,'clim',[-1 1])
title('Original')


subplot(122)
imagesc(squeeze(saturnFilt(:,:,1)))
axis xy, axis square, axis off
set(gca,'clim',[-1 1])
title([ 'Mean-smoothed with k=' num2str(k) ])


%% add large-amplitude spikes
%{
This code loads an image called "saturn.png", converts it to a grayscale image, 
adds Gaussian noise to a random subset of pixels, and applies a median filter to 
the noisy image to reduce the noise. The filtered image is then displayed alongside 
the original and noisy images using subplots.
%}
% we need to load in the image again because I didn't save the original
% before adding Gaussian noise to it (bad programming!).
saturn = double( imread('saturn.png') )./255;

% let's just use the first color dimension for simplicity
saturn = squeeze(saturn(:,:,1));

% add noise to a random subset of pixels
spikelocs = randsample(1:numel(saturn), round(.1*numel(saturn))); % for octave: randsample in the statistics package
%{
This line of code generates a random sample of indices for spiking pixels in the image.
numel(saturn) returns the total number of elements in the image matrix.
round(.1*numel(saturn)) calculates 10% of the total number of elements 
in the matrix and rounds it to the nearest integer.
randsample function is then used to randomly select indices from 1 to numel(saturn) 
without replacement, and the number of indices selected is the 10% of the total number of elements calculated above.

Therefore, spikelocs will contain a random subset of indices for spiking pixels in the image.
%}
noisySaturn = saturn;
noisySaturn(spikelocs) = 123456789;
%{
This line of code sets the values of the elements of noisySaturn at the locations 
specified by spikelocs to the value 123456789. In other words, it adds spikes of 
intensity 123456789 to a random subset of pixels in the noisySaturn image. 
This is done to simulate the presence of random spikes or outliers in the image data.
%}
figure(3), clf, colormap bone

% show the original
subplot(131)
imagesc(saturn)
axis xy, axis square, axis off
set(gca,'clim',[0 1])
title('Original')

% show the noisified version
subplot(132)
imagesc(noisySaturn)
axis xy, axis square, axis off
set(gca,'clim',[0 1])
title('Noisy image')



% k parameter for median filter
k=9;

% notice how I use different variables to avoid overwriting
saturnFilt = noisySaturn;

for i=k+1:size(saturn,1)-k
    for j=k+1:size(saturn,2)-k
        temp = saturn(i-k:i+k,j-k:j+k);
        saturnFilt(i, j) = median(temp(:));
    end % end loop over rows
end % end loop over columns


% and plot the filtered result. Looks pretty good!
subplot(133)
imagesc(saturnFilt)
axis xy, axis square, axis off
set(gca,'clim',[0 1])
title([ 'Median-smoothed with k=' num2str(k) ])

%% Gaussian-based image smoothing

% parameters for Gaussian
gx = -20:20;
gaus2d = zeros(length(gx));

% equal width in both directions. We might call this an isotropic Gaussian.
sx = 5;
sy = 5;

% create the Gaussian point-by-point
for xi=1:length(gx)
    for yi=1:length(gx)
        gaus2d(xi, yi) = exp( -(  (gx(xi)^2)/(2*sx^2)  +  (gx(yi)^2)/(2*sy^2) ));   
        % This line of code defines a 2D Gaussian function that takes two input variables xi and yi, and returns the value of the function at those input values.
    end
end

% plot the 2D gaussian
figure(4), clf
imagesc(gaus2d)
axis image
% hint: try using surf. It looks a bit like a... traffic cone?
%{
imagesc is used to create 2D color-coded images of data. It is often used to display images, 
but can also be used to visualize other types of data. imagesc creates a grid of colored squares 
(pixels) where the color of each square represents the value of the data at that location. 
It can be used to create both continuous and discrete color maps.
%}

% Now smooth the spike-field coherence image.
% (Hint: copy code from previous chapters.)
figure(5), clf
subplot(121)
imagesc(timevec, [], spikeLFP)
set(gca, 'clim', [-1000 1000])

subplot(122)
smo = conv2(spikeLFP, gaus2d,'same');
imagesc(timevec,[],smo./sum(gaus2d(:)))
set(gca,'clim',[-500 500])

%% filter Saturn 
%{
This code performs Fourier analysis on an image of Saturn and applies low-pass and high-pass filters to it. 
Here's a brief overview of what the code does:

1) The image of Saturn is loaded in and its 2D fast Fourier transform (FFT) is computed.
2) A low-pass filter and a high-pass filter are created. The low-pass filter is a rectangular 
function with a width of 2*nPix2use pixels in both the x and y dimensions, centered on 
the middle of the image. The high-pass filter is simply the complement of the low-pass filter.
3) The power spectrum of the image is plotted, with the low-pass and high-pass filters displayed 
in a separate subplot.
4) The image is filtered by multiplying its FFT with the Fourier transform of the filter 
(with the filter centered by applying a shift), and then taking the inverse FFT of the result.
5) The low-pass and high-pass filtered images are displayed in separate subplots.
6) The low-pass filter in this code allows low frequency components to pass through, while blocking 
high frequency components. The high-pass filter does the opposite. The cutoff frequency of the filters 
is controlled by the size of the rectangular function used to create the filter. In this code, the size is set to 2*nPix2use, 
so larger values of nPix2use will result in a lower cutoff frequency.
%}

% load in a fresh copy, just in case
saturn = double( imread('saturn.png') )./255;

% get FFT (remember: 2D!). For simplicity,
% we use only the first color dimension.
saturnX = fft2(squeeze( saturn(:,:,1) ));

% get sizes of image and midpoints
imgdims = size(saturnX);
midX = round(imgdims(2)/2);
midY = round(imgdims(1)/2);

% size of the filter in pixels.
nPix2use = 100;

% create low-pass filter kernel
loPass2d = zeros(imgdims(1:2));
loPass2d(midY-nPix2use:midY+nPix2use,midX-nPix2use:midX+nPix2use) = 1;

% create high-pass filter kernel
hiPass2d = ones(imgdims(1:2));
hiPass2d(midY-nPix2use:midY+nPix2use,midX-nPix2use:midX+nPix2use) = 0;


% and let's see how it looks. First, the (shifted) power spectrum
figure(6), clf
subplot(221)
imagesc(fftshift(log(abs(saturnX))));
axis off, axis square
title('Full power spectrum')


% next, the filter
subplot(222)
imagesc(loPass2d)
axis off, axis square
title('Low-pass filter')

% apply the low-pass filter
filtimg = real(ifft2( saturnX.*fftshift(loPass2d) ));
subplot(223)
imagesc(filtimg)
axis off, axis square
title('Low-pass image')

% and the high-pass filter
filtimg = real(ifft2( saturnX.*fftshift(hiPass2d) ));
subplot(224)
imagesc(filtimg)
axis off, axis square
set(gca,'clim',[-.2 .2]/10)
title('High-pass image')

%% end
