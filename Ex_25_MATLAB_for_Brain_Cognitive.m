%% 25.4.2
%{
Write additional viewing code to show the structural MRI scan in 12
linearly spaced slices between the first and the last, in a 3-by-4 subplot
matrix. Produce three figures corresponding to the three dimensions.
%}
% Load and read in files
cd(['C:\Users\dun81\OneDrive\Documents\Dun II\ET Online\Library Online\Psychology\Computational Neuroscience\' ...
    'MATLAB for Brain and Cognitive Scientists - The MIT Press (2017)\MATLAB_BrainAndCognitiveScientists-main\ch25\'])
[strMRI, pixdim, rotate, dtype] = readnifti('MNI152_T1_1mm.nii');

% create 12 linearly spaced slices on the 3rd dimention of strMRI
linsp = round(linspace(size(strMRI, 3), 1, 12));

figure(1), clf
for ii = 1:length(linsp)
    subplot(3, 4, ii);
    imagesc(squeeze(strMRI(:, :, linsp(ii)))');
    colormap gray
    axis image off
    title(['Slice ', num2str(linsp(ii))]);
end

%% 25.4.4
%{
The structural scan is high resolution. Try down-sampling the image by
creating new images that contain skipped voxels (e.g., take every second
voxel, or every third voxel, and so on). Compute the percentage of
the size of each down-sampled version relative to the original version.
Plot the slice closest to 50% on the z-plane for the different downsampled
images with the proportion of decreased size in the title of the
image.
%}
%{ 
To skip every second voxel, create a loop of linsp that extracts the data
as before (every 12 slices). Then create an if statement within that loop
that says "if at every odd number, process that image, otherwise don't
process that image."
%}
figure(2), clf
for ii = 1:length(linsp)
    for jj = 1:length(linsp)
        if mod(jj, 3) == 1 % check if jj is odd
            jj = ii*jj;
            subplot(3, 4, jj);
            imagesc(squeeze(strMRI(:, :, linsp(jj)))');
            colormap gray
            axis image off
            title(['Slice ', num2str(linsp(jj))]);
        else
            continue
        end
    end
end
%{
The answer shows something simpler. Instead of creating 2 loops, it creates
a single loop that goes from 1:6 instead of 1:12. 
Then similar code to the previous.
Then divide the number of elements (numel function) of the new result to
the original data. Multiply the result by 100 and divide by 100, round them
up and get the result of the percentage. 
%}

%% 25.4.5
%{
Using the fMRI data, compute the mean and standard deviation of the
signal intensity values over time at each voxel. Then show a few slices
from these two maps in a figure. For bonus points, complete this exercise
without using the functions mean and std (or var, that would be
too easy).
%}