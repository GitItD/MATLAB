%% Matlab for Brain and Cognitive Scientists
% Matlab code for Chapter 25
% Mike X Cohen
% 
% This code accompanies the book 
%      "Matlab for Brain and Cognitive Scientists" (MIT Press, 2017).
% Using the code without following the book may lead to confusion, 
% incorrect data analyses, and misinterpretations of results. 
% Mike X Cohen assumes no responsibility for incorrect use of this code. 

%% load and plot structural MRI data
%{
This code is loading and plotting a structural MRI data using the readnifti function. 
It loads the file named "MNI152_T1_1mm.nii" and then displays three slices from 
different planes using the imagesc function.
%}
% load structural mri file
cd(['C:\Users\dun81\OneDrive\Documents\Dun II\ET Online\Library Online\Psychology\Computational Neuroscience\' ...
    'MATLAB for Brain and Cognitive Scientists - The MIT Press (2017)\MATLAB_BrainAndCognitiveScientists-main\ch25\'])
[strMRI pixdim rotate dtype] = readnifti('MNI152_T1_1mm.nii');
%{
This line of code is loading a 3D structural MRI file in NIfTI format with the 
filename 'MNI152_T1_1mm.nii' using the custom function 'readnifti'. 
The function returns several outputs:

'strMRI': a 3D matrix representing the MRI data
'pixdim': a vector containing the voxel dimensions
'rotate': a rotation matrix
'dtype': a data type code for the MRI data.

This code assumes that the 'readnifti' function is available in the current 
MATLAB environment and that the 'MNI152_T1_1mm.nii' file is located in the 
current working directory.
%}
% image difference slices
figure(1), clf
subplot(221)
imagesc(squeeze(strMRI(60, :, :))')
axis image, axis xy

subplot(222)
imagesc(squeeze(strMRI(:, 180, :))')
axis image, axis xy

subplot(223)
imagesc(squeeze(strMRI(:, :, 50))')
axis image, axis xy
%{
The squeeze() function is used to remove any singleton dimensions from the strMRI array, 
so that we can display the 2D slice as an image using the imagesc() function. 
The argument squeeze(strMRI(50,:,:)) takes the 50th slice along the x-axis, and displays 
it as a grayscale image using the imagesc() function. The transpose .' is used to ensure 
that the image is oriented correctly, with the x-axis going from left to right and the y-axis 
going from top to bottom.

The axis image command is used to make sure that the aspect ratio of the image is preserved. 
The axis xy command is used to make sure that the orientation of the y-axis is correct.

Overall, this code is useful for visualizing a single slice of the 3D structural MRI data 
loaded from the file.
%}

colormap gray
set(gcf, 'color', 'k')

% Note: If you look carefully, you can see the axes around the brains. 
% They are dark grey. You can turn them off with
axis off % this command applies only to the last-used axis

%% load in fMRI data
%{
This code is loading fMRI data from multiple files (presumably a sequence of volumes) 
and storing them in a 4D matrix fmridat, where the first three dimensions correspond 
to the x, y, and z voxel coordinates, and the fourth dimension corresponds to time 
(i.e., the volume number in the sequence). It is also extracting the volume number 
from the filename and storing it in the volnums array.
%}
filz = dir('*img');
%{
This line of code uses the dir function to find files in the current directory that have the ".img" extension. 
It is searching for fMRI data files that end in ".img".
%}
for imgi=1:length(filz)
    
    % read one volume
    tempdat = readnifti(filz(imgi).name);
    %{
    The .name property is used to extract the name of the file in the imgith element of the filz struct array. 
The imgi variable is used as an index to access the imgith element in the filz array, and .name is used to 
extract the name of the file stored in that element.

For example, if imgi is 3, filz(imgi).name would give the name of the file stored in the 3rd element 
of the filz struct array.
    
    For examples of how to use this format of coding, see the very bottom
    of this code.
    %}
    % initialize
    if imgi==1 % initialization only during first image!
        fmridat = zeros([size(tempdat), length(filz)]);
        volnums = zeros([1, length(filz)]);
    end
    %{
    zeros([ size(tempdat) length(filz) ]) creates a multi-dimensional array of zeros with the same size as tempdat, 
and with an additional dimension equal to the number of files in filz.

size(tempdat) returns the size of tempdat as a row vector containing the number of rows, columns, and slices of 
the 3D volume, while length(filz) returns the number of files in filz. Therefore, zeros([ size(tempdat) length(filz) ]) 
creates an empty 4D matrix with dimensions [rows, columns, slices, number of files in filz], to which each loaded 
volume is added.

This code creates an fmridat variable and a volnums variable and initializes them to all zeros only 
during the first iteration of the for loop. The if imgi==1 statement ensures that this initialization occurs 
only once, at the beginning of the loop.

The fmridat variable is a 4D matrix with the same size as the tempdat variable (which represents a single 
volume of fMRI data), and with a length in the fourth dimension equal to the number of files that match 
the pattern '*img' in the current directory. The volnums variable is a 1D vector with a length equal to the 
number of files that match the pattern '*img', and all elements are initialized to zero.
    %}
    % put volume into larger matrix
    fmridat(:, :, :, imgi) = tempdat;
    
    
    % and get file name
    uscore = strfind(filz(imgi).name, '_');
    dotloc = strfind(filz(imgi).name, '.');
    volnums(imgi) = sscanf(filz(imgi).name(uscore+1:dotloc-1), '%g');
    %{
    This code extracts the volume number from the filename of each fMRI image file in the directory.

The function strfind is used to find the index of the first occurrence of the character '_'  in the filename, 
which separates the volume number from the rest of the filename. The variable uscore stores the index of '_'.

The function strfind is also used to find the index of the first occurrence of the character '.' in the filename, 
which separates the file extension from the rest of the filename. The variable dotloc stores the index of '.'.

The function sscanf is then used to extract the volume number from the portion of the filename between '_' and '.'. 
The format specifier '%g' is used to indicate that the number should be interpreted as a floating-point value. 
The extracted volume number is then stored in the volnums array.
    %}
end

%% plot a few images and time courses

figure(2), clf, set(gcf,'color','w')

% image one slice (try some others!)
subplot(221)
imagesc(squeeze(fmridat(20, :, :, 1))')
%{
The code imagesc(squeeze(fmridat(20,:,:,1))') is trying to display the 20th slice of the first volume of the fMRI data fmridat. 
The squeeze function is used to remove any singleton dimensions in the 4D matrix, so that fmridat has dimensions of [X Y Z T]. 
The imagesc function is then used to display the 2D image data of the 20th slice of the first volume. 
The ' at the end is used to transpose the matrix before displaying it.
%}
axis xy, axis image, axis off
set(gca,'clim',[0 12000])
title('Volume 20')

% same slice but later in time
subplot(222)
imagesc(squeeze(fmridat(20, :, :, end-10))')
axis xy, axis image, axis off
set(gca,'clim', [0, 12000])
title([ 'Volume ' num2str(size(fmridat, 4)-10) ])
%{
size(fmridat, 4) is a MATLAB function that returns the size of the fourth dimension of the fmridat array. 
In this code, it is being used to get the number of volumes in the fmridat array, assuming that the volumes 
are stored in the fourth dimension of the array.

(size(fmridat, 4)-10) is a mathematical expression that subtracts 10 from the number of volumes in the fmridat array.
Subtracting 10 from this value gives the index of the volume to be displayed in the title of the plot.
%}
colormap gray

% and a time course
subplot(212)
plot(squeeze(fmridat(20, 40, 40, :)))
title(sprintf('Time course of voxel I, J, K = %g, %g, %g', 20, 40, 40))
%{
This code sets the title of the current figure to a string that includes the voxel coordinates of the time course being plotted. 
The sprintf function is used to create the string using a format specifier, which includes placeholders for the voxel coordinates. 
The %g placeholder specifies that the corresponding argument (20, 40, and 40) should be formatted as a floating-point number. 
The placeholders are replaced by the values in the order they are listed after the format specifier.

sprintf is a built-in MATLAB function that allows you to create a formatted string from input arguments. 
It takes a string as the first argument, and then one or more additional arguments that specify the values to be 
inserted into the string. The values can be of any MATLAB data type, and the format of the output can be controlled 
using various format specifiers.

The three %g are placeholders for the three floating-point numbers that are being passed as arguments to the sprintf function.
%}
xlabel('Time (volume)')
ylabel('Signal intensity (a.u.)')

%% experiment design

onsets = (6:12:84)+12;

timeline = zeros(length(filz), 1);
for i=0:5
    timeline(onsets+i) = 1;
end
%{
This code creates an array onsets with the starting points of each block of events, spaced 12 volumes apart. 
It then initializes an array timeline with zeros for each volume in the dataset. The for loop iterates through each 
block of events (there are 6 blocks) and sets the corresponding volumes in timeline to 1, representing the onset of the event block. 
The onsets+i index selects each volume within the current event block, which is spaced 1 volume apart.
%}
%% numerator of t-test

m0 = mean(fmridat(:, :, :, timeline==0), 4);
m1 = mean(fmridat(:, :, :, timeline==1), 4);

numerator = m1-m0;
%{
mean(fmridat(:, :, :, timeline==0), 4) calculates the mean of the fMRI data over the fourth dimension (time) only 
for those timepoints where timeline is 0. Similarly, mean(fmridat(:, :, :, timeline==1), 4) calculates the mean of the 
fMRI data over the fourth dimension only for those timepoints where timeline is 1.

Then, numerator calculates the difference between the mean fMRI data for the two conditions (m1-m0), voxel-wise 
for all the three spatial dimensions, resulting in a 3D image. This difference image can be used to identify brain regions 
that show differential activity during the presence and absence of the stimulus.
%}
%% denominator of t-test

v0 = var(fmridat(:, :, :, timeline==0), [], 4);
v1 = var(fmridat(:, :, :, timeline==1), [], 4);

denominator = sqrt( (v0/sum(timeline==0)) + (v1/sum(timeline==1)) );
%{
This code computes the denominator of the t-statistic used for voxel-wise statistical inference in functional MRI data analysis.

The variance (v0, v1) of each voxel in the two groups (here defined by the timeline variable) is computed using the var function 
applied along the fourth dimension (i.e., the time dimension) of the fmridat array. The optional second input argument, [], specifies 
that the variance should be computed using the default normalization by N-1, where N is the sample size. The resulting output has the 
same dimensions as the input fmridat, except for the fourth dimension which is collapsed to a single value.

The denominator is then computed as the square root of a weighted sum of the variances. The weights are the reciprocal of the sample
 sizes in each group (i.e., the number of time points in each condition, computed as the sum of the timeline variable for each group). 
The division by the sum of the timeline variable in each group normalizes the weights so that they sum to 1. 
This weighted sum is equivalent to a pooled estimate of the standard deviation across the two groups.
%}
%% t-statistic

% map of t-statistic values
tmap = numerator ./ denominator;

% fairly arbitrarily threshold at t<2.5
tthresh = tmap;
tthresh(tthresh<2.5) = 0;

% let's have a look at a thresholded slice
figure(3), clf
imagesc(squeeze(tthresh(:, :, 35))')

%{
This code generates a map of t-statistic values, by dividing the difference between the mean activation levels 
of the two conditions by the standard deviation of the data across the two conditions. 
Then, the code sets a threshold value for t-statistics and sets all values below that threshold to zero. 
Finally, it displays a single slice of the thresholded t-map.

Here's a breakdown of the code:

tmap = numerator ./ denominator;: This line calculates the t-statistic map by dividing the 
difference in means by the standard deviation of the data across the two conditions.
tthresh = tmap; tthresh(tthresh<2.5) = 0;: This line sets a threshold value of 2.5 for 
the t-statistic map and sets all values below this threshold to zero.
figure(3), clf; imagesc(squeeze(tthresh(:, :, 35))'): This line displays a single slice of the 
thresholded t-statistic map, in this case the 35th slice of the data along the third dimension (i.e. the axial plane). 
figure(3) opens a new figure window, and clf clears any previous contents from the figure. 
imagesc is used to display the thresholded t-map, and squeeze is used to remove any singleton dimensions from the data. 
The transpose operator ' is used to display the image in the correct orientation.
%}

%% threshold 1-voxel results

% slightly better threshold based on contiguously 
% significant regions. bwconncomp is in the image
% processing toolbox.
islands = bwconncomp(tthresh);
% find the number of pixels in each 
islandsizes = cellfun(@length, islands.PixelIdxList);

for ii=1:islands.NumObjects
    if islandsizes(ii)<3
        tthresh(islands.PixelIdxList{ii}) = 0;
    end
end
%{
This code uses the MATLAB function bwconncomp from the Image Processing Toolbox to identify 
contiguous regions (islands) of significant t-values in the tthresh map generated in the previous code block.

The function bwconncomp takes a binary image as input and returns a structure islands containing information 
about the connected components of the binary image. Each connected component is represented by a list of 
pixel indices in the PixelIdxList field of the islands structure.

The code then calculates the size of each connected component using the cellfun function, which applies 
the function @length to each element of the cell array islands.PixelIdxList.

Finally, the code iterates through each connected component and sets the t-values to 0 if the component has 
fewer than 3 pixels. This is a simple thresholding step that eliminates very small islands that are likely to be noise.
%}
%% now plot

figure(4), clf

% pick a slice to show (fix to be axial)
slice2plot = 37;

% Color discritization. 64 is the Matlab standard.
% Try other numbers to see what happens to the plot. 
% Try, for example, 5 vs. 500
cdiscr = 64;

% extract slice and normalize image to [0 64]
img2plot = squeeze(m1(:,:,slice2plot))';
img2plot = img2plot - min(img2plot(:));
img2plot = cdiscr * img2plot./max(img2plot(:));


% same for stats map
stat2plot = squeeze(tthresh(:,:,slice2plot))';
stat2plot(stat2plot==0) = NaN;
% this time, normalize image to [64 128] because it will 
% be plotted on top of the structural map.
stat2plot = stat2plot - min(stat2plot(:));
stat2plot = cdiscr+cdiscr*.5 + cdiscr*.5 * stat2plot./max(stat2plot(:));

% uncomment the next line to reproduce the white voxels in the book figure
% stat2plot(stat2plot>1) = max(stat2plot(:));

% use pcolor function, which is similar to surf.
% Try replacing pcolor with surf and turning on rotate3d
h(1) = pcolor( img2plot );
hold on
h(2) = pcolor( stat2plot );
axis xy, axis image, axis off

colormap([ gray(cdiscr); hot(cdiscr) ]);
set(h,'linestyle','none')
%{
This code generates a plot of an MRI image and a statistical map overlayed on top of the image.

The first section of the code (lines 5-14) sets up the variables needed to calculate the statistical map. 
m0 and m1 are the means of the fMRI data for the baseline and task conditions, respectively, and numerator 
is the difference between the two means. v0 and v1 are the variances of the fMRI data for the baseline and task 
conditions, respectively, and denominator is the square root of the sum of the variances divided by the number 
of data points for each condition.

The next section of the code (lines 16-22) generates a thresholded version of the statistical map. 
Voxels with a t-value less than 2.5 are set to 0.

The final section of the code (lines 24-45) creates the plot. The img2plot variable is the MRI image for the chosen slice. 
It is normalized to have values between 0 and 64. The stat2plot variable is the thresholded statistical map for the chosen slice. 
It is normalized to have values between 64 and 128 and is plotted on top of the MRI image using the pcolor function. 
The resulting plot has a colorbar that goes from black to white to yellow, with black representing the lowest intensity 
and yellow representing the highest intensity.
%}
%% end.

%{

% create a struct array with two fields
person(1).name = 'Alice';
person(1).age = 25;
person(2).name = 'Bob';
person(2).age = 30;

% access the fields using the dot notation
disp(person(1).name); % output: Alice
disp(person(2).age); % output: 30

% create a struct array with three fields
fruit(1).name = 'apple';
fruit(1).color = 'red';
fruit(2).name = 'banana';
fruit(2).color = 'yellow';
fruit(3).name = 'orange';
fruit(3).color = 'orange';

% loop through the struct array and display the names and colors
for i = 1:length(fruit)
    disp(fruit(i).name);
    disp(fruit(i).color);
end

%}