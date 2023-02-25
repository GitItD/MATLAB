%% 16.7.1
%{
In this exercise, you will explore the effects of correlation strength and
sample size on simulated correlated data. Write a triple-loop to create
two signals of different lengths and correlation strengths and 100 iterations
(trials) of each length-strength pair. For each simulated correlation,
compute the squared distance between the true correlation that
you specified and the estimated correlation coefficient. Plot the trialaveraged
results in two 2D images (one image for average correlation
coefficients, and one image for the average distance to the true correlation.
How do you interpret the results?is it better to have a strong
effect or more data?
%}
ii = 10;
jj = 50;
nn = 100;
data = randn(nn, 2);
dataa = randn(nn, 2);
r = rand()*0.9 + 0.1;
% Write a triple-loop to create 2 signals of different lengths and
% correlation strengths:
for ii = 1:ii
    data(:, ii) = data(:,1)*r(:,ii) + data(:,2)*sqrt(1 - r(:,ii)^2);
    data2(:,ii) = corr(data(:,ii));
    for jj = 1:jj
        dataa(:, jj) = dataa(:,1)*r(:,jj) + dataa(:,2)*sqrt(1 - r(:,jj)^2);
        dataa2(jj) = bsxfun(@divide, mean(dataa(:,jj)), std(dataa(:,jj)));
        for nn = 1:nn
            c(nn) = bsxfun(@minus, dataa2(:, nn), data2(:, nn));
        end
    end
end
% 100 trials of each length-strength pair. For each simulated correlation,
% compute square distance between true correlation and estimated
% correlation coefficient

% Plot trialaveraged results in two 2D images

%Answer by ChatGPT:
% Define correlation strengths and signal lengths
strengths = [0.2 0.5 0.8];
lengths = [50 100 150];

% Initialize arrays to store results
avg_corrs = zeros(length(strengths), length(lengths));
avg_dists = zeros(length(strengths), length(lengths));

% Loop over strength-length pairs
for i = 1:length(strengths)
    for j = 1:length(lengths)
        % Initialize arrays to store results for this pair
        corrs = zeros(1, 100);
        dists = zeros(1, 100);
        
        % Simulate 100 trials for this pair
        for k = 1:100
            % Generate two signals of different lengths and correlation strengths
            x = randn(lengths(j), 1);
            y = randn(lengths(j)+10, 1);
            y = y(1:lengths(j)) + strengths(i) * x;

            % Compute the correlation coefficient
            corrs(k) = corr(x, y);

            % Compute the squared distance between the true correlation and the estimated correlation
            dists(k) = (strengths(i) - corrs(k))^2;
        end
        
        % Compute the trial-averaged correlation coefficient and distance to true correlation
        avg_corrs(i, j) = mean(corrs);
        avg_dists(i, j) = mean(dists);
    end
end

% Plot the trial-averaged results
figure(1), clf
imagesc(avg_corrs)
title('Average Correlation Coefficients')
xlabel('Signal Length')
ylabel('Correlation Strength')
colorbar

figure(2), clf
imagesc(avg_dists)
title('Average Distance to True Correlation')
xlabel('Signal Length')
ylabel('Correlation Strength')
colorbar

% Answer by author:
ns = round(linspace(10, 200, 60));
cs = linspace(.1, .9, 40);
nreps = 100;

corrmat = zeros(length(ns), length(cs), nreps, 2);

for ni = 1:length(ns)
    for ci = 1:length(cs)
        for ti = 1:nreps
            
            x = randn(ns(ni), 2);
            x(:,2) = x(:,1)*cs(ci) + x(:,2)*sqrt(1-cs(ci)^2);
            
            x = bsxfun(@minus, x, mean(x, 1));
            x = bsxfun(@rdivide, x, std(x, [], 1));
            tmp = x'*x/(ns(ni)-1);
            corrmat(ni, ci, ti, 1) = tmp(1, 2);
            
            corrmat(ni, ci, ti, 2) = (tmp(1, 2) - cs(ci))^2;
        end
    end
end

figure(3), clf
subplot(211)
contourf(cs, ns, squeeze(mean(corrmat(:,:,:,1),3)), 40, 'linecolor', 'none')
set(gca, 'clim', [0 1])

subplot(212)
contourf(cs, ns, (squeeze(mean(corrmat(:,:,:,2), 3))), 80, 'linecolor', 'none')
set(gca, 'clim', [0, 0.025])
xlabel('Correlation strength'), ylabel('N')

%% 16.7.9
%{
Compute the channels-to-channels covariance matrix that was used in
figure 16.4 in two different ways. First, average all trials together before
computing the covariance matrix. Second, compute the covariance
matrix for each trial separately and then average together 99 covariance
matrices (there are 99 trials). Plot the two results next to each
other, using the same color scaling. Are there striking differences, and
how would you explain them?
%}

% Get data and compress dimensions for calculation readiness
load(['C:\Users\dun81\OneDrive\Documents\Dun II\ET Online\Library Online\Psychology\' ...
    'Computational Neuroscience\MATLAB for Brain and Cognitive Scientists - The MIT Press (2017)\MATLAB_BrainAndCognitiveScientists-main\sampleEEGdata.mat'])
data = squeeze(EEG.data(:,:,1));

% Average all trials together
data_avg = mean(data(:));

% Compute covariance matrix separately
% Method 1
cov1 = data*data'/(EEG.pnts - 1);   
%Method 2
data = bsxfun(@minus, data, mean(data, 2));
cov2 = data*data'/(EEG.pnts - 1);

% Average together 99 covariance
cov1avg = mean(cov1(:));
cov2avg = mean(cov2(:));

% Plot the results next to each other
figure(4), clf
set(gca, 'clim', [0 200])
colormap hot
colorbar
title('Covariance Comparisons')

subplot(121)
topoplotIndie(diag(cov1), EEG.chanlocs);

subplot(122)
topoplotIndie(diag(cov2), EEG.chanlocs);

% Author's answer:
load(['C:\Users\dun81\OneDrive\Documents\Dun II\ET Online\Library Online\Psychology\' ...
    'Computational Neuroscience\MATLAB for Brain and Cognitive Scientists - The MIT Press (2017)\MATLAB_BrainAndCognitiveScientists-main\sampleEEGdata.mat'])

% covariance of average
data2 = squeeze(mean(EEG.data,3));
data2 = bsxfun(@minus,data2,mean(data2,2));
covAve = data2*data2'/(EEG.pnts-1);

% covariance of single-trials
covTrl = zeros(EEG.nbchan);
for ti=1:EEG.trials
    data2 = squeeze(EEG.data(:,:,ti));
    data2 = bsxfun(@minus,data2,mean(data2,2));
    covTrl = covTrl + data2*data2'/(EEG.pnts-1);
end

covTrl = covTrl./EEG.trials;

figure(5), clf
subplot(121)
imagesc(covAve), set(gca,'clim',[-10 10]), colorbar
title('Covariance of average')

subplot(122)
imagesc(covTrl), set(gca,'clim',[-120 120]), colorbar
title('Average of single-trial covariances')

% Are the answers the same? Let ChatGPT decide.
%{
No, the two codes are not equivalent. They are calculating and displaying different covariance matrices.

The first code computes two covariance matrices (cov1 and cov2) using different methods 
and then averages them separately, whereas the second code computes two different types of 
covariance matrices (covAve and covTrl) and displays them in two separate subplots of a figure.

covAve is the covariance matrix of the average of all trials, while covTrl is the average of the 
covariance matrices of each individual trial.
%}

% To make the first code equivalent to the second one, you can modify it as follows:
% Load the data and extract the relevant data
load(['C:\Users\dun81\OneDrive\Documents\Dun II\ET Online\Library Online\Psychology\' ...
    'Computational Neuroscience\MATLAB for Brain and Cognitive Scientists - The MIT Press (2017)\MATLAB_BrainAndCognitiveScientists-main\sampleEEGdata.mat'])
data = squeeze(EEG.data);

% Compute the covariance of the average
data_avg = mean(data, 3);
data_avg = bsxfun(@minus, data_avg, mean(data_avg, 2));
covAve = data_avg*data_avg'/(EEG.pnts - 1);

% Compute the average of the single-trial covariances
covTrl = zeros(EEG.nbchan);
for ti=1:EEG.trials
    data_trl = squeeze(EEG.data(:,:,ti));
    data_trl = bsxfun(@minus, data_trl, mean(data_trl, 2));
    covTrl = covTrl + data_trl*data_trl'/(EEG.pnts-1);
end
covTrl = covTrl./EEG.trials;

% Plot the results next to each other
figure(6), clf
subplot(121)
imagesc(covAve), set(gca,'clim',[-10 10]), colorbar
title('Covariance of average')

subplot(122)
imagesc(covTrl), set(gca,'clim',[-120 120]), colorbar
title('Average of single-trial covariances')
