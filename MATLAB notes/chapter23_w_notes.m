%% Matlab for Brain and Cognitive Scientists
% Matlab code for Chapter 23
% Mike X Cohen
% 
% This code accompanies the book 
%      "Matlab for Brain and Cognitive Scientists" (MIT Press, 2017).
% Using the code without following the book may lead to confusion, 
% incorrect data analyses, and misinterpretations of results. 
% Mike X Cohen assumes no responsibility for incorrect use of this code. 

%% import data

% import data from previous chapter
cd(['C:\Users\dun81\OneDrive\Documents\Dun II\ET Online\Library Online\Psychology\Computational Neuroscience\' ...
    'MATLAB for Brain and Cognitive Scientists - The MIT Press (2017)\MATLAB_BrainAndCognitiveScientists-main\ch22\'])
fid = fopen('../ch22/03-04-11-aa_sig1_spikes.dat');         % Open the file and assign it to fid.

spikenum = 1;
trialnum = 0;
spikesparse = zeros(100, 2);                   
spikesfull  = zeros(1, 1001);

while ~feof(fid)            % feof = test for end of file is true (1) or false (0)... 
    
    % get a line of code
    d=fgetl(fid);           % fgetl = read line from file, and returns it as a string. Each call to 'fgetl' reads next line and stores it into 'd'.
    
    % if it's an empty line, skip to the next line
    if isempty(d), continue; end
    
    % if this is an inter-trial comment, increment trial number 
    % and then skip forward
    if d(1)=='#'    % This code block checks whether the current line of data (d) starts with a # character. This is commonly used to indicate comments or meta-information in data files.
        if strcmpi(d(1:5),'# int') 
            %{ 
strcmpi is a MATLAB function that compares two strings ignoring the case of the characters. It returns 1 if the two strings are the same, ignoring case, and 0 otherwise.
strcmpi(d(1:5),'# int') = compares the first 5 characters of the string d with the string '# int'.
            %}
            trialnum = trialnum + 1;
        end
        continue;
        
        
    % otherwise, it's real data; grab the spike time 
    % and update the matrices
    else
        spiketime = ceil(sscanf(d,'%g')/1000);
        %{
        sscanf is a function in MATLAB that reads formatted data from a string. In this case, sscanf(d,'%g') is 
reading the numerical values from the string d, which is a line of text from a file. The %g format specifier 
tells sscanf to read the numerical values in decimal.

For example, if d is the string "10.5 20.3 30.1", then sscanf(d,'%g') will return the vector [10.5, 20.3, 30.1]. 
The values are automatically converted to floating point numbers.
        %}
        spikesparse(spikenum,:) = [trialnum, spiketime];
        spikesfull(trialnum, spiketime) = 1;
        spikenum = spikenum+1;
    end
end
fclose(fid);
%{
This block of code updates two matrices spikesparse and spikesfull based 
on the spike times in the imported data file.

spikesparse is a sparse matrix with two columns where each row corresponds to a 
spike time in a specific trial. The first column contains the trial number and the second 
column contains the spike time in milliseconds.

spikesfull is a full matrix where each row corresponds to a trial and each column corresponds 
to a time bin (in this case, 1 millisecond). If a spike occurred in a given time bin for a given trial, 
the corresponding entry in spikesfull is set to 1, otherwise it is 0.
The spikesfull matrix is initialized to all zeros using the line spikesfull = zeros(1,1001); before the data is read in. 
Then, in the else branch of the if statement where d(1)~='#', spikesfull(trialnum, spiketime) = 1; assigns a value 
of 1 at the appropriate location in the matrix where the spike occurred. Therefore, the elements of spikesfull 
will be either 0 or 1 depending on whether or not a spike occurred at that time and trial.

spikenum is a counter variable that keeps track of the number of spikes that have been processed 
so far, while trialnum keeps track of the trial number that is associated with the spikes being processed.

spikesfull(trialnum, spiketime) = 1;: This line updates the spikesfull matrix by setting the value of the entry 
corresponding to the current trial and spike time to 1.

spikenum = spikenum+1;: This line increments the spikenum counter so that the next spike will be 
added to the next row of the spikesparse matrix.
%}

figure(1), clf
plot(spikesparse(:,2), spikesparse(:,1), 'm.')
xlabel('Time (ms)'), ylabel('Trials')

%% spike timing with full matrices
%{
This code computes the spike rate over time for a given window around each spike. 
It first defines a window size of 50 ms, and initializes an empty vector spikerhyth to hold the spike rates.
%}
win = 50; % in ms and also in indices (only in this dataset!)... time window of spike count
spikerhyth = zeros(1, win*2+1); % vector of zeros... used to sum up the spikes in each time window. Equal to win*2+1 b/c we want the vector to center around 0.
n = 0; % counter variable used to calculate average spike number in each window.
%{
The reason for initializing the spikerhyth variable with zeros is to create an array of the 
same size and data type that will hold the spike counts.

In this code, spikerhyth is used to accumulate the number of spikes that occurred at each 
time bin in a sliding window. If you didn't initialize it with zeros before adding the spike counts, 
the array would be empty, and you wouldn't be able to add the spike counts to it.

By initializing the array with zeros, you are effectively creating a "canvas" of the correct size and 
shape that you can then draw the spikes onto. In other words, it's like preparing a blank piece of 
paper before you start drawing a picture on it.
%}
for triali = 1:trialnum
    
    spikes = find(spikesfull(triali,:)); % Find all timepoints where there are spikes by finding indices of nonzero elements.
    % need an extra line of code here...
    
    for si = 1:length(spikes)
        spikerhyth = spikerhyth + spikesfull(triali, spikes(si) - win:spikes(si) + win); % Each time window of zeros is added to the spike points
        n = n + 1; % counter variable from earlier
    end % end spike loop
end % end trial loop
%{
This loop iterates over each trial in the spikesfull matrix (which has dimensions trialnum by 1001, 
where each row represents a different trial and each column represents a different millisecond).

For each trial, it finds all the time points (in milliseconds) where there were spikes by finding 
the indices of nonzero elements using the find function. Then, for each spike time point, it adds 
the spike window (of length win milliseconds) to the spikerhyth vector using indexing, which 
represents the summed spike raster within a time window around each spike. It also increments 
n by 1 for each spike, which will be used later to calculate the average spikerhyth.

Note that there is a missing line of code between the spikes and si loops. This line should 
increment triali by 1. Without this line, the loop will always execute on triali = 1 and will not iterate 
over any other trials.
%}
% divide by N to finalize average
spikerhyth = spikerhyth./n;
%{
This line of code is dividing each value in the spikerhyth vector by the total number of spikes n in all trials. 
This normalization step finalizes the average spike rate at each time bin of the spikerhyth vector.
%}
figure(2), clf
plot(-win:win, spikerhyth, 'rs-')    % Plot the average spike rate at each time bin. The -win:win represents the beginning and end of time window of the plot, and spikerhyth is the each timepoint where the spike data exists.
set(gca, 'xlim', [-10 10])

%% spike timing with sparse matrices
%{
This section of code is computing spike timing histograms using a sparse matrix format. 
The sparse matrix format saves memory by only storing the non-zero values of the matrix, 
so it is more efficient when most of the matrix values are zero.
%}
win = 50; % in ms and also in indices
spikerhyth = zeros(1, win*2+1); % vector of zeros... used to sum up the spikes in each time window. Equal to win*2+1 b/c we want the vector to center around 0.
n = 0; % counter variable used to calculate average spike number in each window.

for triali=1:max(spikesparse(:,1)) 
    %{ 
    max(spikesparse(:,1)) returns the maximum value in the first column of the spikesparse matrix. 
Since spikesparse contains spike times and their corresponding trial numbers, the maximum value 
in the first column represents the highest trial number in the matrix. In other words, this expression 
returns the total number of trials in the experiment.
    %}
    spikes = spikesparse(spikesparse(:,1) == triali, 2);  
    %{
    selects the 2nd column of spikesparse (which contains spike times) for all rows where the trial 
number in the 1st column matches the current trial number triali. These spike times are stored in the spikes variable.
    %}
    
    for si=1:length(spikes)
        tempspikes = spikes-spikes(si);
        tempspikes(tempspikes < -win | tempspikes > win) = [];
        
        spikerhyth(tempspikes + win + 1) = spikerhyth(tempspikes + win + 1) + 1;
        
        n = n+1;
        
        %{
        This code calculates a spike time histogram using sparse matrix. 
For each trial, it first extracts the spike times from the spikesparse matrix where 
the first column corresponds to the trial number and the second column corresponds to 
the spike time. Then, for each spike in the extracted spike times, it calculates the time differences 
between the current spike time and all other spikes in the same trial using the tempspikes = spikes - spikes(si) line.

Next, it filters the tempspikes such that it only retains those values that fall within a time window of [-win, win].
 This is achieved by setting to [] any tempspikes that are less than -win or greater than win.

Finally, it adds one to the element in spikerhyth corresponding to each non-zero tempspikes. 
The n variable keeps track of the number of spikes used to generate the histogram.
        %}
        
    end % end spike loop
end % end trial loop
spikerhyth = spikerhyth./n;     % normalize spikerhyth


% plot on top of the previous result to see how they compare
figure(2), hold on
plot(-win:win, spikerhyth, 'ko-')
set(gca,'xlim', [-10 10])
legend({'From full matrix'; 'From sparse matrix'})

%% using autocorrelation in the frequency domain
%{
This code computes the autocorrelation of the spike train using the power spectrum of the Fourier transform. 
The process starts by initializing an array spikerhyth with zeros of length win*2+1. Then, for each trial, it first 
subtracts the mean of the trial from the data. It then takes the FFT of the mean-subtracted data, and computes 
the power spectrum by squaring the absolute value of the FFT. An inverse FFT is applied to the power spectrum 
to obtain the autocorrelation function. To account for edge effects, the relevant parts of the autocorrelation function 
are selected by indexing the last win+1 points of the end and first win+1 points. These parts are then summed into 
spikerhyth. After iterating through all the trials, spikerhyth is divided by the total number of trials. Finally, the result is 
plotted as a red diamond marker in Figure 3, with a size of 10. The x-axis limits are set to [-10 10].
%}
spikerhyth = zeros(1, win*2+1);

for triali=1:trialnum
    
    % mean subtract data from this trial
    meansub = bsxfun(@minus, spikesfull(triali,:), mean(spikesfull(triali,:)));
    
    % take FFT and iFFT of power spectrum
    temp = fft(meansub);
    temp = ifft( abs(temp).^2 );
    
    % chop out parts we want
    temp = temp([end-win+1:end 1:win+1]);
    
    % sum autocorrelation
    spikerhyth = spikerhyth + temp./max(temp);
end

% again, divide by N
spikerhyth = spikerhyth./triali;


% and plot
figure(3), clf
plot(-win:win,spikerhyth,'rd-','markersize',10)
set(gca,'xlim',[-10 10])

%% the previous cell can also be done more efficiently
%{
This code uses a more efficient way to calculate the same result as the previous cell.

The long line of code inside the spikerhyth variable is doing the following:

bsxfun(@minus, spikesfull, mean(spikesfull,2)): it subtracts the mean of each row of spikesfull from each element in that row.
fft(...,[],2): it computes the Fast Fourier Transform (FFT) of each row of the previous result.
abs(...).^2: it takes the absolute value of each element of the previous result and squares it.
ifft(...,[],2): it computes the inverse FFT of each row of the previous result.
mean(...,1): it calculates the mean of all the rows of the previous result, which is the same as adding up all the rows and then dividing by the number of rows.
The resulting spikerhyth vector is the same as in the previous cell.
%}
% Using one really long line of code! So I guess it depends on your definition of "efficient."
spikerhyth = mean( ifft( abs(fft( bsxfun(@minus,spikesfull,mean(spikesfull,2)),[],2 )).^2,[],2 ) ,1) / trialnum;

figure(3), hold on
plot(-win:win, spikerhyth([end-win+1:end 1:win+1]), 'm')
set(gca, 'xlim', [-10 10])

%% cross-neuron spike timing correlation
% In this code block, we are analyzing the spike timing correlation between two different neuron clusters.

load times_090425blk10_ch115.mat

% spike timing with sparse matrices
win = 50; % in ms
spikerhyth = zeros(1, win*2+1);
n = 0;

spikesClus1 = cluster_class(cluster_class(:,1) == 1,2);
spikesClus2 = cluster_class(cluster_class(:,1) == 2,2);
%{ 
Extract variable from cluster_class(:,1) that is equal to
cluster_class(1,2) and cluster_class(2, 2), then put it into
spikesClus1 and spikesClus2 respectively.
%}

for si = 1:length(spikesClus1) % create a loop with length of the data
    tempspikes = round(spikesClus2-spikesClus1(si)); 
    tempspikes(tempspikes<-win | tempspikes>win) = [];
    %{
    This line calculates the time differences between spikes in spikesClus1 and spikesClus2. 
The differences are rounded to the nearest integer and stored in the tempspikes variable.
    %}
    spikerhyth(tempspikes+win+1) = spikerhyth(tempspikes+win+1)+1;
    %{
    This line increments the values in the spikerhyth array based on the elements in tempspikes. 
Specifically, it increments the values at indices tempspikes+win+1 by 1.
    %}
    n = n+1;
end % end spike loop
spikerhyth = spikerhyth./n; % normalize spikerhyth

figure(4), clf
plot(-win:win,spikerhyth,'kd-','markersize',10)
set(gca,'xlim',[-10 10])
xlabel('Relative spike times (ms)')
ylabel('Probability')

%% hc-5 for spike-lfp coherence
%{
This code loads data from the hc-5 dataset, which contains spike and local field potential (LFP) recordings 
from rat hippocampus. It then performs some basic analysis and plotting on the data
%}

% data are taken from:
% https://crcns.org/data-sets/hc/hc-5/about-hc-5/
load spikefieldData

wherespikes = find(spikeTimes); % finds the indices in spikeTimes where there are spikes
win = 400; %  sets the size of the analysis window to 200 indices (not milliseconds). The window will be centered around each spike.
wherespikes(wherespikes<win | wherespikes>length(lfp)-win) = []; % removes any spikes that are too close to the beginning or end of the lfp recording to allow for a full analysis window.
spikeLFP = zeros(length(wherespikes), win*2+1); % initializes the spikeLFP variable as a matrix with one row for each spike and a column for each index in the analysis window.

for si=1:length(wherespikes)
    spikeLFP(si,:) = lfp(wherespikes(si)-win:wherespikes(si)+win); % This extracts the LFP data for the current spike and stores it in the corresponding row of spikeLFP.
end

figure(5), clf
timevec = (-win:win)/(1000/srate);
plot(timevec, mean(spikeLFP))
hold on
plot([0 0], get(gca,'ylim'),'k')
set(gca, 'xlim', timevec([1 end]))
xlabel('Peri-spike time (ms)'), ylabel('Voltage (mV)')
%{
The code then plots the mean LFP waveform across all spikes in spikeLFP, along with a vertical line at time 0 (the spike time). 
The x-axis is labeled "Peri-spike time (ms)" even though the units are actually indices.
%}

figure(6), clf
subplot(121)
h = plot(timevec, spikeLFP(1:100:end,:));
hold on
plot(timevec, mean(spikeLFP), 'k', 'linew', 4)
set(gca, 'xlim', timevec([1 end]))
set(h, 'color', [.7 .7 .7])
xlabel('Peri-spike time (ms)'), ylabel('Trials')
%{
The code then creates a subplot with two panels (subplot(121) and subplot(122)). 
In the left panel, the code plots the LFP waveform for every 100th spike in spikeLFP, 
along with the mean LFP waveform across all spikes.
%}

subplot(122)
imagesc(timevec, [], spikeLFP)
set(gca, 'clim', [-500 500], 'xlim', timevec([1 end]))
xlabel('Peri-spike time (ms)'), ylabel('Trials')
%{
The right panel uses the imagesc function to create a heatmap of the LFP data for all spikes in spikeLFP. 
The x-axis is labeled "Peri-spike time (ms)" even though the units are actually indices, and the y-axis is labeled "Trials". 
The color scale is set to range from -500 to 500.
%}

%% same but filtered
%{
Here, a complex Morlet wavelet is being created in order to filter the LFP signal. 
The wavelet is defined by a center frequency (frex) and a number of cycles (nCycl). 
The wavelet is defined in the time domain and is convolved with the LFP signal in the frequency domain.
%}
frex = 8; % hz

wavetime = -2:1/srate:2;
halfwave = floor(length(wavetime)/2)+1;


nData = length(lfp);
nWave = length(wavetime);
nConv = nData + nWave - 1;

% create wavelet
nCycl = 6;
gausS = nCycl / (2*pi*frex);

wavelet  = exp(2*1i*pi*frex*wavetime + (-wavetime.^2)/(2*gausS^2) );
waveletX = fft(wavelet, nConv);
waveletX = waveletX./max(waveletX); % normalize

%% convolution
%{
This code computes the spike-LFP coherence using a wavelet-based approach.
Overall, the code provides a way to visualize the relationship between spikes 
and the LFP signal at a specific frequency.

The first line computes the Fourier transform of the LFP signal lfp using the 
fft function with a specified length nConv. This could be used to smooth the 
data or to remove unwanted frequency components.

The second line convolves the Fourier transform of the LFP data with the 
wavelet transform waveletX, which is presumably a complex-valued vector. 
This produces a complex-valued output as representing the amplitude and 
phase of the LFP signal at each scale (or frequency) of the wavelet.
%}
dataX = fft(lfp, nConv);
as = ifft( dataX .* waveletX );

% cut off edges
%{
The third line removes a half wave length of the wavelet transform from 
each end of the as vector. This is likely done to avoid edge effects caused 
by the convolution of the LFP signal with the wavelet transform.
%}
as = as(halfwave-1:end-halfwave);

% now
%{
In this line of code, angle is the MATLAB function that computes the angle 
(in radians) of a complex number. as(logical(spikeTimes)) selects the entries 
of as where spikeTimes is true, and these entries are treated as complex numbers 
with their real part being the original value of as and their imaginary part being zero.
%}
sfc = abs(mean(exp(1i*angle(as(logical(spikeTimes))))));

wherespikes = find(spikeTimes); % find places where spikes have occured
win = 200; % indices, not ms!
wherespikes(wherespikes<win | wherespikes>length(as)-win) = [];
%{
A window size of 200 is defined and any spikes that occur within the first win or last win 
time points are excluded. The length of as is used as the upper limit for the last win exclusion.
%}
spikeLFP = zeros(length(wherespikes),win*2+1); % Matrix to be filled; 1 row of length 'win*2+1' for each spike.

for si=1:length(wherespikes)
    spikeLFP(si,:) = real(as(wherespikes(si)-win:wherespikes(si)+win)); % assign spikeLFP with real values of continuous wavelet transform 'as'
end

figure(7), clf
timevec = (-win:win)/(1000/srate);
plot(timevec, mean(spikeLFP))
set(gca, 'xlim', timevec([1 end]))
%{
This code plots the average local field potential (LFP) around spikes. 
It first creates a time vector timevec representing the time relative to each 
spike, using the win variable which determines the size of the time window 
in milliseconds. It then plots the mean of the spikeLFP matrix, which contains 
the LFP data surrounding each spike. The set function sets the x-axis limits to the 
first and last values of timevec.

Note that this code assumes that the sampling rate (srate) is in Hz.
%}

figure(8), clf
subplot(121)
h = plot(timevec, spikeLFP(1:100:end,:));
hold on
plot(timevec, mean(spikeLFP), 'k', 'linew', 4)
set(gca, 'xlim', timevec([1 end]))
%{
set(gca, 'xlim', timevec([1 end])) sets the limits of the x-axis of the current axis (gca) to the 
minimum and maximum values of the timevec vector, which is the time axis for the plotted data. 
In other words, it sets the x-axis limits to the first and last time points of the plotted data.
%}
set(h, 'color', [.7 .7 .7])

subplot(122)
imagesc(timevec, [], spikeLFP)
set(gca, 'clim', [-200 200])
%{
In this case, imagesc is being used to visualize the spike-LFP data matrix as an image, 
where each row represents a spike-triggered LFP waveform. The function imagesc 
takes as input the x-axis values (in this case, timevec), the y-axis values (empty brackets, 
meaning the row indices will be automatically determined), and the matrix to be visualized (spikeLFP). 
The matrix values are mapped to a color scale, where the minimum value is mapped to the lowest color 
and the maximum value is mapped to the highest color. By default, the color scale ranges from blue to red, 
with white representing the midpoint.
%}

%% now smooth the image

gx = -10:10; 

gaus2d = zeros(length(gx)); % create window of zeros

s = 10;

for xi=1:length(gx) % count to length of window
    for yi=1:length(gx) % another count equal to xi
        gaus2d(xi,yi) = exp( -( gx(xi)^2+gx(yi)^2 ) / (2*s^2) ); % create 2D Gaussian filter for smoothing image. Both xi and yi increased until gx are multiplied
    end
end


figure(9), clf

subplot(211)
imagesc(gaus2d)
axis image

subplot(212)
smo = conv2(spikeLFP, gaus2d, 'same');
%{
The conv2 function in MATLAB is used for 2-dimensional convolution. 
It computes the full convolution of matrices A and B, i.e., the output is the 
size of the sum of the sizes of the input matrices minus one in each dimension. 
The 'same' option in the conv2 function returns the central part of the convolution 
that is the same size as matrix A.

In the example code provided, the spikeLFP matrix is convolved with gaus2d using 
the 'same' option to obtain the smoothed image. This convolution is performed element-wise 
on each 2D submatrix in spikeLFP. The resulting matrix smo has the same size as the spikeLFP matrix.
%}
imagesc(timevec, [], smo./sum(gaus2d(:)))
set(gca, 'clim', [-100 100])
%{
The code imagesc(timevec, [], smo./sum(gaus2d(:))) generates an image where each row 
corresponds to a different spike, and the columns represent the smoothed LFP at different 
time points around the spike. The third argument smo./sum(gaus2d(:)) scales the smoothed 
LFP values by the sum of the Gaussian filter values so that the resulting image represents an 
estimate of the LFP surrounding the spikes with the Gaussian smoothed filter.

The timevec argument sets the horizontal axis to represent time in seconds. The empty brackets in imagesc(timevec,[],...) are used to specify that the vertical axis should be auto-scaled to the size of the data being plotted.
%}
%%