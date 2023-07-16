%% Matlab for Brain and Cognitive Scientists
% Matlab code for Chapter 22
% Mike X Cohen
% 
% This code accompanies the book 
%      "Matlab for Brain and Cognitive Scientists" (MIT Press, 2017).
% Using the code without following the book may lead to confusion, 
% incorrect data analyses, and misinterpretations of results. 
% Mike X Cohen assumes no responsibility for incorrect use of this code. 

%% analog vs. binary representation of a spike

% use data from chapter 23
cd(['C:\Users\dun81\OneDrive\Documents\Dun II\ET Online\Library Online\Psychology\Computational Neuroscience\' ...
    'MATLAB for Brain and Cognitive Scientists - The MIT Press (2017)\MATLAB_BrainAndCognitiveScientists-main\ch22\'])
load(['03-04-11-aa_sig1_spikes.dat'])
%{
The focus of the analyses in this chapter is on data recorded during tasks
with repeated discrete trials; for example, when data were recorded while
experimental animals viewed images, experienced somatosensory stimulation,
or initiated motor responses. The data weÅfll use in this chapter
were downloaded from crcns.org (https://crcns.org/data-sets/motor-cortex/
alm-1/). In the first data set (Li et al. 2015), a mouse was instructed to stick
his tongue out repeatedly.

This code loads data from a file named 03-04-11-aa_sig1_spikes.dat and plots 
it in two subplots to show the difference between analog and binary representation of a spike.
%}
figure(1), clf

% analog
%{
The first subplot shows the analog representation of the voltage data of the 10th spike trial. 
The x-axis represents time in arbitrary units (a.u.) and the y-axis represents voltage in millivolts (mV).
%}
subplot(211)
plot(spikes(10,:),'k','linew',2)
set(gca,'xlim',[0 13])
xlabel('Time (a.u.)')
ylabel('Voltage (mV)')

% binary simplification
%{
The second subplot shows the binary representation of the same data. 
The x-axis and y-axis labels are the same as in the first subplot, but the y-axis 
shows only two levels - 0 (no spike) and 1 (spike). The black lines indicate the 
time of the spike, which is simplified to a binary representation.
%}
subplot(212)
set(gca, 'xlim', [0 13], 'ylim', [-.1 1.2], 'ytick', [0 1], 'yticklabel', {'no';'yes'})
hold on
plot([1 12], [0 0],'k','linew',2)
plot([5 5], [0 1],'k','linew',2)
xlabel('Time (a.u.)')
ylabel('Spike?')

%% import data

%{
Our first goal will be to import the data into MATLAB from the data files
provided online. We will put the data into two matrices, one full and one
sparse.

Viewing this file with a text editor reveals
that the file contains some meta-data followed by spike times, and then
more meta-data, and so on. Each section of meta-data is indicated by hash
marks (#) and separates different trials of data. We donÅft know a priori how
much data there will be, so a while-loop seems appropriate.
%}
cd(['C:\Users\dun81\OneDrive\Documents\Dun II\ET Online\Library Online\Psychology\Computational Neuroscience\' ...
    'MATLAB for Brain and Cognitive Scientists - The MIT Press (2017)\MATLAB_BrainAndCognitiveScientists-main\ch22\'])
fid = fopen('03-04-11-aa_sig1_spikes.dat');

%{
The spike times will be stored in two matrices: APsS is for Ågaction potentials
sparseÅh and APsM is for Ågaction potentials matrix.Åh APsS will be a sparse
matrix that stores the trial number and time stamp of each spike; APsM will
be a full trials-by-time matrix in which each element of the matrix is 0 (no
spike) or 1 (spike).
%}
spikenum = 1;
trialnum = 0;
APsS = zeros(100,2);
APsM = zeros(1,1001);

%{
We are ready to read in the data from the file. The while-loop will
continue until the end of the file; in other words, while ~feof(fid). The
code inside the while-loop is separated into four sections. The first two
sections are fairly straightforward: read in a line of code, and skip to the
next line of code if the current line is empty.
%}
while ~feof(fid)    % feof is a MATLAB function that returns true when the end of a file has been reached, and false otherwise.
    
    % get a line of code
    d=fgetl(fid);
    %{
    fgetl is a MATLAB function that reads a line of text from an open file identified by a file identifier (a numeric value returned 
by the fopen function). It returns the line of text as a string, excluding the line terminator. If the end of the file is reached, fgetl returns -1.
    %}
    
    % if it's an empty line, skip to the next line
    if isempty(d), continue; end
    
    % if this is an inter-trial comment, increment trial number and then skip forward
    %{
    The goal of the next section of the code is to skip through the trial information
section. If the line starts with a hash, we can skip forward. But we
want special treatment for the last line of the comments, because that is our
cue to increment the trial count variable.
    
    If you use a text editor to look into the data file, you will see that each
trial contains information about the stimulus intensity. The full line of text
is Ågintensity (dB): 97.3333Åh but there is no other line that starts with Ågint,Åh
which means we need to match only the first few characters. If there is no
match, the if-then statement is false, trialnum is not updated, and the
script skips forward to the next line of the file.
    %}
    if d(1) == '#'
        if strcmpi(d(1:5),'# int')      % strcmpi is a MATLAB function that compares two strings for equality ignoring the case of the characters. It returns 1 if the two strings are equal, and 0 otherwise. For example, strcmpi('hello', 'HELLO') returns 1.
            trialnum = trialnum + 1;
        end
        continue;
        %{
        This block of code checks whether the current line of text starts with the character #. 
If it does, it checks whether the first five characters of the line are equal to # int, 
which indicates the start of a new trial. If so, it increments the trialnum variable by 1. 
Finally, the continue statement skips to the next iteration of the while-loop, effectively skipping 
over the current line of text. This block of code is used to skip through the meta-data sections of 
the file and only process the spike time data.
        %}
        
    % otherwise, it's real data; grab the spike time and update the matrices
    %{
    The fourth section of the while-loop extracts the spike time and enters it
into the matrices. Notice the difference between how the data are entered
in the sparse versus the full action potential matrices.
    
    This section of code is executed when the current line of code is not an 
empty line or a comment line. It extracts the spike time from the line of code using 
sscanf, which reads formatted data from a string. Here, %g is used as the format 
specifier to read a floating-point number, and ceil is applied to convert 
the spike time to an integer value.
    %}
    else
        spiketime = ceil(sscanf(d, '%g')/1000);
        %{
        sscanf is a MATLAB function that reads formatted input from a string. It reads 
data from the string according to the format specified by the format specifier and 
returns the data in one or more output variables. In the code snippet above, sscanf(d, '%g') 
reads the spike time from the current line of the file fid and returns it as a double-precision number. 
The %g format specifier reads a floating-point number, including scientific notation, in decimal or 
hexadecimal notation. The ceil function is then used to round up the spike time to the nearest millisecond, 
which is the time resolution used in the data.
        
        ceil is a built-in MATLAB function that returns the smallest integer greater than or equal to the input 
argument. For example, ceil(3.2) returns 4, ceil(2.8) returns 3, and ceil(-2.3) returns -2.
        %}
        APsS(spikenum,:) = [trialnum spiketime];
        APsM(trialnum,spiketime) = 1;
        spikenum = spikenum + 1;
    end
end
fclose(fid);

npnts = size(APsM, 2);
timevec = 1:npnts;
%{
The first line sets the variable npnts to the number of columns (i.e., time points) in the matrix APsM.
The second line creates a vector called timevec which contains the values 1 through npnts. 
This is a time vector that can be used as the x-axis when plotting the spike data.
%}

% plot the data as an image
figure(2), clf
subplot(311)
imagesc(bsxfun(@times, APsS',[10 1]')), colormap gray
xlabel('Spike count')
set(gca,'ytick',1:2,'yticklabel',{'Trial';'Time'})

% and a different kind of image
subplot(312)
plot(APsS(:,2),APsS(:,1),'k.','markersize',3)
set(gca,'xlim',timevec([1 end]),'ylim',[1 trialnum])
xlabel('Time (ms)'), ylabel('Trial')

% ...the same kind of image but using a different format
subplot(313)
imagesc(timevec, [], 2-APsM)
xlabel('Time (ms)'), ylabel('Trial')

%% spike rate (spikes/second)
%{
So far weÅfve computed the number of spikes in each trial between 0 and
1,000 milliseconds. To make this finding more comparable across different
studies, stimuli, and so forth, the number should be converted to units of
spikes per second (often written sp/s).

That is, although 10 Hz and 10 sp/s both mean that the neuron emitted
10 action potentials within a period of 1 second, Åg10 HzÅh implies that
there was one spike roughly each 100 milliseconds, whereas Åg10 sp/sÅh simply
indicates that a total of 10 spikes were observed somewhere in the time
frame of 1 second.

To convert from spike count to spikes per second, simply divide by the
total amount of time in seconds
%}
tidx = dsearchn(timevec', [0 1000]');

dt = timevec(tidx(2)-tidx(1)) / 1000;
spikeRates = sum(APsM(:,tidx(1):tidx(2)),2) / dt;

figure(3), clf
plot(spikeRates)
xlabel('Trial'), ylabel('spike rate (sp/s)')
set(gca, 'xlim', [0 size(APsM,1)+1])

%% PSTH

%{
Extracting the average spike count within a time window is useful but provides
coarse temporal information. You might want to know how the average
spike rate changes over time. The standard way to examine the time
course of spiking activity is to make the time windows shorter, and then
have many windows over time. The results produce a histogram-like plot of
averaged spike rates over time.

This code generates two plots, one for the peri-stimulus time histogram (PSTH) 
at the original temporal resolution of the data and one for the PSTH downsampled 
to a smaller number of bins.
%}

figure(4), clf

% PSTH at the temporal resolution of the data
%{
In the first subplot (211), the mean spike rate of all trials is calculated for each 
time point, and plotted against the time vector. The time vector is set as the x-axis 
label and the y-axis label represents the spike rate.
%}
subplot(211)
dt = mean(diff(timevec)) / 1000;
plot(timevec,mean(APsM)./dt,'k','linew',2)
set(gca,'xlim',timevec([1 end]))
xlabel('Time (ms)'), ylabel('Spike rate')

% downsampling to a smaller number of bins
%{
In the second subplot (212), the data is downsampled by summing the spikes over 
bins of 5 ms width. The time vector and the corresponding spike rate for each bin are 
calculated and plotted against each other. The time vector is set as the x-axis label and 
the y-axis label represents the spike rate. The mean spike rate is calculated over 
all trials in each bin.
%}
subplot(212)
bins = ceil(timevec/5);
[spPerbin,timevecBin] = deal(zeros(1, max(bins)));
for i=1:length(spPerbin)
    spPerbin(i) = mean( sum(APsM(:,bins==i), 2) ,1);
    timevecBin(i) = mean(timevec(bins==i));
end

dt = mean(diff(timevecBin))/1000;
plot(timevecBin,spPerbin./dt,'k','linew',2)
set(gca,'xlim',timevecBin([1 end]))
xlabel('Time (ms)'), ylabel('Spike rate')
%{
The dt variable is calculated as the mean time difference between consecutive time 
bins in the timevecBin variable, divided by 1000 to convert from milliseconds to seconds. 
This is then used to normalize the spike rate calculated in each bin (stored in the spPerbin variable), 
so that the units of the y-axis are in spikes per second (sp/s).

The downsampling of the data into a smaller number of bins (in this case, using 5 ms bins) is done 
to reduce the amount of noise in the data and to make it easier to visualize patterns in the spike 
rate over time.
%}

%% end.
