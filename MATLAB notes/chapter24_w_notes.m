%% Matlab for Brain and Cognitive Scientists
% Matlab code for Chapter 24
% Mike X Cohen
% 
% This code accompanies the book 
%      "Matlab for Brain and Cognitive Scientists" (MIT Press, 2017).
% Using the code without following the book may lead to confusion, 
% incorrect data analyses, and misinterpretations of results. 
% Mike X Cohen assumes no responsibility for incorrect use of this code. 

%% simple 2D categorization
%{
Spike sorting notes are found in "Mathematics of Action Potential and Spike
Trains".

First, we will test to see whether if we can seperate signals of spikes
from different groups of neurons by amplitude and widths of the spikes:
We will define spike amplitude as the distance between the peak of the
spike waveform and the preceding dip. And we will define spike width as
the distance in time between the pre-peak dip and the post-peak dip.


%}
% use data from previous chapter
cd(['C:\Users\dun81\OneDrive\Documents\Dun II\ET Online\Library Online\Psychology\Computational Neuroscience\' ...
    'MATLAB for Brain and Cognitive Scientists - The MIT Press (2017)\MATLAB_BrainAndCognitiveScientists-main\ch23\'])
load ../ch23/times_090425blk10_ch115.mat

spikefeat = zeros(size(spikes, 1), 2);

for spikei=1:size(spikes,1)
    
    % find peak
    [~,peakidx] = max(spikes(spikei,:));
    
    % find minimum before peak
    [~,min1idx] = min(spikes(spikei,1:peakidx));
    
    % find minimum after peak
    [~,min2idx] = min(spikes(spikei, peakidx:end));
    min2idx = min2idx+peakidx-1;
    
    % get premin-peak difference
    spikefeat(spikei,1) = diff( spikes(spikei,[min1idx peakidx]) );
    
    % get min2min
    spikefeat(spikei,2) = min2idx-min1idx;
end

%{
Now letfs see if these two features can differentiate different neurons. We
will make a plot of one feature versus the other feature (spike amplitude vs.
spike width), with each dot reflecting one action potential. If all action
potentials come from one neuron, there will be little spatial structure in the
plot. However, if the action potentials recorded by this electrode come from
two different neurons, there will be separable clouds of dots, with each
cloud corresponding to each neuron.
%}

% plot one feature by the other
figure(1), clf
h = plot(spikefeat(:,1),spikefeat(:,2),'ko','markerfacecolor','k');
set(gca,'ylim',[2 12],'xlim',[-5 300])
set(h,'markeredgecolor',[.7 .7 .7])
xlabel('Spike amplitude')
ylabel('Spike width')

%{
Discussion about the results:
The resulting plot do not allow us to distinguish between two neuron
populations and which signal was attributed to each type. There are two
possibilities for why this is:
1. All the spikes in this dataset comes from one group of neurons. 
2. Spike amplitude and widths are terrible parameters for distinguishing
between two populations of neuron groups. 
%}

%% now with PCA
%{
Amplitude and width alone are insufficient to isolate two neuron datasets.
We can apply PCA to waveforms and inspect whether the projection of each
spike onto principal component (PC) space can differentiate two neuron
sets.

PCA involves eigendecomposition of a covariance matrix.
%}
% subtract mean
spikes = bsxfun(@minus, spikes, mean(spikes, 2));

% covariance matrix
spikecov = spikes'*spikes / (size(spikes, 1) - 1);     % Here if we put A'*A, we are taking the dot product, which results in square matrix of each column of A column x A column, rather than square matrix of size A row x A row (which is the case for A*A')
%{
(spikes, 1) is passing the argument 1 as the second input to the mean function, 
which specifies the dimension along which the mean is calculated. In this case, 
1 indicates that the mean should be calculated across the rows (i.e., for each time 
point across all spikes).

(size(spikes, 1) - 1) calculates the degrees of freedom for the sample covariance matrix calculation. 
The sample covariance matrix is an estimate of the true covariance matrix of the population, calculated 
using a finite number of samples. The degrees of freedom represent the number of independent observations
in the sample that are used to estimate the covariance. In this case, subtracting 1 from the number of samples 
ensures that the sample covariance matrix is an unbiased estimator of the true covariance matrix.
%}
[eigvects, eigvals] = eig(spikecov);

figure(2), clf
subplot(221)
plot(mean(spikes, 1))
set(gca, 'xlim', [.75 12.25], 'ylim',[-45 85])
xlabel('Time (a.u.)')
title('All spikes')

subplot(223)
plot(eigvects(:,end-1:end),'linew',2)
legend({'PC1'; 'PC2'})
set(gca,'xlim',[.75 12.25],'ylim',[-.7 .7])
title('Component time courses')
xlabel('Time (a.u.)')

subplot(222)
imagesc(spikecov)
set(gca,'clim', [-3000 3000])
xlabel('Time point'), ylabel('Time point')
title('Covariance matrix')

%% now with ICA

r = rank(spikecov);

weights = jader(spikes', r);     % jader is apparently supposed to be a function, but causes an error. This is because I don't have the toolbox installed. I've learned how to install new tools, the notes can be found in OneNote-->E-Learn-->MATLAB Scripts and Functions
icas = weights*spikes';

figure(3), clf
plot(icas(1,:),icas(2,:),'k.','markersize',5)
xlabel('IC #1'), ylabel('IC #2')
set(gca,'xlim',[-7 0],'ylim',[-6 6])

% for Octave, kmeans is in the statistics package
ICclustidx = kmeans(icas(1:2,:)',2);

figure(4), clf
plot(icas(1,ICclustidx==1),icas(2,ICclustidx==1),'r.','markersize',3), hold on
plot(icas(1,ICclustidx==2),icas(2,ICclustidx==2),'b.','markersize',3)
xlabel('IC #1'), ylabel('IC #2')
legend({'cluster 1';'cluster 2'})

%% project each spike onto first and second components
%{
Generating PCA plots. The codes below clearly shows us two distinct groups
of neuron signals.
%}
comp1 = spikes*eigvects(:,end);
comp2 = spikes*eigvects(:,end-1);


figure(5)
plot(comp1, comp2,'k.','markersize',5)
xlabel('PC #1'), ylabel('PC #2')
% The axis limits were adjusted for the book, although it cuts off a few
% spikes. same for the next figure, and for the ICA-based analyses.
set(gca,'xlim', [-280 -0], 'ylim', [-120 120])
%

% separate into two groups (you'll learn about kmeans in chapter 31).
PCclustidx = kmeans([comp1 comp2], 2);

figure(6), clf
h1 = plot(comp1(PCclustidx==1), comp2(PCclustidx==1),'r^','markersize',4); hold on
h2 = plot(comp1(PCclustidx==2), comp2(PCclustidx==2),'bo','markersize',4);
xlabel('PC #1'), ylabel('PC #2')
set(gca,'xlim',[-280 -0],'ylim',[-120 120])
legend({'cluster 1';'cluster 2'})

set(h1,'color','k')
set(h2,'color',ones(1,3)*.5)

%% end.
