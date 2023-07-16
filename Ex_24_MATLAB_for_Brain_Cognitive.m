%% 24.5.4
%{
Are PCA-based and ICA-based clusters the same? What is the proportion
of spikes that overlap? Remember that labeling groups as Åg1Åh or
Åg2Åh is arbitrary; youÅfll need to devise a method to see if they match
even if PCA-derived unit Åg1Åh corresponds to ICA-derived unit Åg2.Åh For
the two clusters, plot the average waveform shape for the overlapping
PCA and ICA clusters, and, in a separate plot, show the difference
between PCA and ICA component time courses. Next, compare these
to the clusters that were already formed in the downloaded data set
(look for the variable called cluster_class).
%}
cd(['C:\Users\dun81\OneDrive\Documents\Dun II\ET Online\Library Online\Psychology\Computational Neuroscience\' ...
    'MATLAB for Brain and Cognitive Scientists - The MIT Press (2017)\MATLAB_BrainAndCognitiveScientists-main\ch23\'])
load ../ch23/times_090425blk10_ch115.mat
% Create PCA clusters
spikes = bsxfun(@minus, spikes, mean(spikes, 2));
spikecov = spikes'*spikes / (size(spikes, 1) -1); 
[eigvects, eigvals] = eig(spikecov);

% Create ICA clusters
r = rank(spikecov);
weights = jader(spikes', r);  
icas = weights*spikes';

% Measure the proportion of spikes that overlap
comp1 = spikes*eigvects(:,end);
comp2 = (icas(end-1,:))';

% Plot the average waveform shape for overlapping PCA and ICA clusters
figure(1), clf
plot(comp1, comp2,'k.','markersize',5)
xlabel('PC #1'), ylabel('IC #2')

% Show the difference between PCA and ICA component time courses
PCclustidx = kmeans([comp1 comp2], 2);

figure(2), clf
h1 = plot(comp1(PCclustidx==1), comp2(PCclustidx==1),'r^','markersize',4); hold on
h2 = plot(comp1(PCclustidx==2), comp2(PCclustidx==2),'bo','markersize',4);
xlabel('PC #1'), ylabel('IC #2')
set(gca,'xlim',[-280 -0],'ylim',[-120 120])
legend({'cluster 1';'cluster 2'})

set(h1,'color','k')
set(h2,'color',ones(1,3)*.5)

% Compare the results to clusters that were already formed class_class
figure(3), clf
plot(mean(cluster_class, 2), 'k.','markersize',5)
xlabel('PC #1'), ylabel('IC #2')
%{
figure(3), clf
for ii = 1:size(cluster_class, 2)
    if cluster_class(1, ii) == 1
        plot((cluster_class(:, ii)), 'r^', 'markersize', 4); hold on
    else
        plot((cluster_class(:, ii)), 'bo', 'markersize', 4);
    end
end

xlabel('PC #1'), ylabel('PC #2')
set(gca,'xlim',[-280 -0],'ylim',[-120 120])
legend({'cluster 1';'cluster 2'})

set(h1,'color','k')
set(h2,'color',ones(1,3)*.5)
%}
%% 24.5.5
%{
Neurons have a refractory period and cannot emit an action potential
too closely after the preceding action potential. If a spike cluster is
really from a single neuron, there should be no double-spikes within a
few milliseconds of each other. In large data sets, there will be a small
number of such fast double-spikes, which may reflect noise or simultaneous
spikes from multiple neighboring neurons. In the two clusters
here, compute the interspike-interval histograms and show in a plot.
Write code to count the number of spikes that occur within 1 millisecond,
2 milliseconds, and so on, up to 20 milliseconds, after each spike,
separately for the two clusters.
%}

% Assign each spike to PCA or ICA cluster
pca_cluster_spikes = comp1;
ica_cluster_spikes = comp2;

% Compute interspike intervals
pca_cluster_isi = diff(pca_cluster_spikes);
ica_cluster_isi = diff(ica_cluster_spikes);

% Bin intervals into histograms
isi_bins = 0:1:20; % define bin edges
pca_cluster_hist = histogram(pca_cluster_isi, isi_bins);
ica_cluster_hist = histogram(ica_cluster_isi, isi_bins);

%{ 
Plot histograms
figure;
subplot(2,1,1);
bar(pca_cluster_hist.BinEdges(1:end-1), pca_cluster_hist.Values);
xlabel('Interspike interval (ms)');
ylabel('Count');
title('PCA cluster ISI histogram');

subplot(2,1,2);
bar(ica_cluster_hist.BinEdges(1:end-1), ica_cluster_hist.Values);
xlabel('Interspike interval (ms)');
ylabel('Count');
title('ICA cluster ISI histogram');

% Count num spikes occur within 1 ms - 20 ms w/ each other for cluster 1:

% Count num spikes occur within 1 ms - 20 ms w/ each other for cluster 2:
%}

%% 24.5.7
%{
Now that you have two spike clusters, test whether those units interact.
Identify action potentials in cluster no. 2 that spike either 3 milliseconds
before or 3 milliseconds after spikes in cluster no. 1. Plot their
waveform shapes. For comparison, also plot the average waveform
shapes of randomly selected spikes.
%}