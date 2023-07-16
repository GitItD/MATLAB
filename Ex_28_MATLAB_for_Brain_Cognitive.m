%% 28.6.12
%{
Run the reaction time EEG correlations again, but apply permutation
testing to transform the coefficients into statistical z-scores. Consult
chapter 15 if you need a refresher on permutation testing for statistical
inference. What do you shuffle at each iteration during permutation
testing? Generate three plots: one showing the raw coefficients from
the least-squares fit, one showing the statistical z-values, and one in
which all z-values with an absolute value less than 1.96 are turned to
zero (z = 1.96 corresponds to p < 0.05, two-tailed and uncorrected for
multiple comparisons).
%}

cd(['C:\Users\dun81\OneDrive\Documents\Dun II\ET Online\Library Online\Psychology\' ...
    'Computational Neuroscience\MATLAB for Brain and Cognitive Scientists - The MIT Press (2017)\MATLAB_BrainAndCognitiveScientists-main'])
load sampleEEGdata.mat

rts = zeros(size(EEG.epoch));

% loop over trials
for i=1:EEG.trials    
    % find the index corresponding to time=0, i.e., trial onset
    [~, zeroloc] = min(abs( cell2mat(EEG.epoch(i).eventlatency) ));    
    % reaction time is the event after the trial onset
    rts(i) = EEG.epoch(i).eventlatency{zeroloc + 1};
end
% create design matrix
A = [ones(EEG.trials, 1), rts'];
x = A(:,2);
y = 2*x + randn(size(x));
yHat = x.*A(:,1) + y.*A(:,2);
% create regression line
b = regress(y, [ones(size(x)), x]);
% setup blank canvas
zscores = zeros(size(b));
y_perm = zeros(size(y));
% fill blank canvas
for i = 1:EEG.trials
y_perm = y(randperm(length(y)));
b_perm = regress(y_perm, [ones(size(x)), x]);
zscores = zscores + abs(b_perm - b);
end

mu_null = mean(zscores);
sigma_null = std(zscores);

% compute the z-score of the original coefficients
z = abs(b - mean(b)) ./ std(b);

% plot the null distribution and the z-score of the original coefficients
figure(1)
subplot(121)
plot(1:length(x), y*500, 'o', 'markerface', 'b', 'markersize', 5)
hold on
plot(1:length(x), yHat)
title('Least-squares Fit')
subplot(122)
histogram(zscores);
hold on;
plot([mu_null mu_null], [0 max(histcounts(zscores))], 'r', 'LineWidth', 2);
plot(z, 0, 'ko', 'MarkerSize', 10, 'LineWidth', 2);
xlabel('Absolute difference in coefficients');
ylabel('Frequency');
legend({'Null distribution', 'Mean of null distribution', 'Z-score of original coefficients'});
% End
%% 28.6.12 ChatGPT
cd(['C:\Users\dun81\OneDrive\Documents\Dun II\ET Online\Library Online\Psychology\' ...
    'Computational Neuroscience\MATLAB for Brain and Cognitive Scientists - The MIT Press (2017)\MATLAB_BrainAndCognitiveScientists-main'])
load sampleEEGdata.mat

rts = zeros(size(EEG.epoch));

% loop over trials
for i=1:EEG.trials    
    % find the index corresponding to time=0, i.e., trial onset
    [~, zeroloc] = min(abs( cell2mat(EEG.epoch(i).eventlatency) ));    
    % reaction time is the event after the trial onset
    rts(i) = EEG.epoch(i).eventlatency{zeroloc + 1};
end

% create design matrix
A = [ones(EEG.trials, 1), rts'];
x = A(:,2);
y = 2*x + randn(size(x));
yHat = x.*A(:,1) + y.*A(:,2);

% create regression line
b = regress(y, [ones(size(x)), x]);

% setup blank canvas
zscores = NaN(size(b));
y_perm = NaN(size(y));

% fill blank canvas
parfor i = 1:EEG.trials
    y_perm = y(randperm(length(y)));
    b_perm = regress(y_perm, [ones(size(x)), x]);
    zscores(i) = abs(b_perm - b);
end

mu_null = nanmean(zscores);
sigma_null = nanstd(zscores);

% compute the z-score of the original coefficients
z = abs(b - mean(b)) ./ std(b);

% plot the null distribution and the z-score of the original coefficients
figure(1)
subplot(121)
plot(1:length(x), y*500, 'o', 'markerface', 'b', 'markersize', 5)
hold on
plot(1:length(x), yHat)
title('Least-squares Fit')

subplot(122)
histogram(zscores);
hold on;
plot([mu_null mu_null], [0 max(histcounts(zscores))], 'r', 'LineWidth', 2);
plot(z, 0, 'ko', 'MarkerSize', 10, 'LineWidth', 2);
xlabel('Absolute difference in coefficients');
ylabel('Frequency');
legend({'Null distribution', 'Mean of null distribution', 'Z-score of original coefficients'});

% set threshold for significant z-scores
threshold = 1.96;

% convert z-scores to zero if their absolute value is less than the threshold
z(abs(z) < threshold) = 0;

% plot the z-scores with thresholding
figure(2)
bar(z)
title('Statistical z-scores')
xlabel('Coefficients')
ylabel('Z-score')

% perform permutation testing and convert the resulting z-scores to zero if their absolute value is less than the threshold
zscores_perm = zeros(size(b));
for i = 1:1000 % run 1000 iterations
    y_perm = y(randperm(length(y)));
    b_perm = regress(y_perm, [ones(size(x)), x]);
    zscores_perm = zscores_perm + abs(b_perm - b - mean(b_perm - b)) ./ std(b_perm - b); % compute z-scores and add them to the sum
end
zscores_perm = zscores_perm / 1000; % divide the sum by the number of iterations to obtain the mean z-score
zscores_perm(abs(zscores_perm) < threshold) = 0; % convert z-scores to zero if their absolute value is less than the threshold

% plot the z-scores from permutation testing with thresholding
figure(3)
bar(zscores_perm)
title('Statistical z-scores from permutation testing')
xlabel('Coefficients')
ylabel('Z-score')

% plot the z-scores from permutation testing with thresholding and zero out insignificant values
figure(4)
z_thresholded = zscores_perm;
z_thresholded(abs(z_thresholded) < threshold) = 0;
bar(z_thresholded)
title('Statistical z-scores from permutation testing with thresholding')
xlabel('Coefficients')
ylabel('Z-score')
% end

%% 28.6.12 Mike X Cohen

%{
This code is an example of time-frequency analysis applied to EEG data, and then correlation analysis is 
performed between time-frequency representation and reaction times. Here is an overview of what this code is doing:

1. The EEG data is loaded from a file called `sampleEEGdata.mat`.
2. The reaction times are computed for each trial based on the onset of the trial and the timing of the relevant event.
3. Wavelets are created in the frequency domain using the `freqrange` and `numfrex` parameters. These wavelets are used to 
decompose the EEG signal into time-frequency representations.
4. The time-frequency representations are computed for each trial by convolving the EEG data with the wavelets for each frequency.
5. The correlation between the time-frequency representation and reaction times is computed using a least-squares regression approach. 
The result is a correlation matrix for each frequency and time point.
6. The significance of the correlation is evaluated using a permutation test that generates fake data based on randomization of reaction times, 
and a Z-score map is created to show the correlation between time-frequency representation and reaction times. The resulting Z-score map is 
displayed as a plot of frequency versus time.
%}

%% 28.6.14
%{
The reaction time EEG correlations involved a design matrix with two
variables (intercept and reaction time). Redo the analysis using a third
independent variable corresponding to beta-band power. The value of
the third column on each trial should be power from 14 to 20 Hz and
from 600 to 200 milliseconds from that trial. How do the results
look? Of course, the correlation in the gseedh window should not be
interpreted. This is one method of analyzing brain functional interactions
(in this specific case, the results suggest that trials with stronger
prestimulus beta-band power also had stronger poststimulus thetaband
power).
%}
%% 28.6.16
%{
(This question is related to chapter 25.) Hemodynamic activity does
not suddenly increase and decrease, like the boxcar-shaped regressor
we used in the t-test. FMRI analyses involve generating a canonical
hemodynamic response and convolving the design matrix (the independent
variable) with that canonical hemodynamic response (the kernel).
Search the Internet to find a function called spm_hrf.m and use it
to generate a canonical hemodynamic response function. Then convolve
the regressor (variable timeline in the code in chapter 25) with
the hemodynamic response function kernel using techniques you
learned in chapter 12. Then create a design matrix that includes an
intercept and the regressor, and fit the model to each voxelfs time
course. How do these results compare with the t-test we implemented
in chapter 25?
%}