clear
clc

%% 15.8
%{
The MATLAB function to compute a percentile score from a distribution
is called prctile. This function is in the Statistics and Machine
Learning toolbox. If you have this toolbox, use the help file for that
function to figure out how to compute the 98% percentile of a distribution
of 10,000 random integers uniformly distributed between
14 and 2,345 (you will also need to generate this distribution).
Next, devise a solution to this problem that is not dependent on the
toolbox.
%}

%% My answer:

% Generate data:
n = (14:2345);
nperm = 10000;
rn = randi([14, 2345], [1, 10000]);

% Generate mean and std:
rnavg = mean(rn);
rnstd = std(rn);

%Calculated SE:
se = rnstd/sqrt(nperm);
se = se*2.33;
% Create normal distribution
norn = normrnd(rnavg, rnstd);
cminus = rnavg - se
cplus = rnavg + se

figure(1), clf
hold on
histogram(rn, 100)

%% ChatGPT's answer:

% Generate 10,000 random integers uniformly distributed between 14 and 2,345
rng(1); % seed the random number generator for reproducibility
dist = randi([14, 2345], 1, 10000);

% Sort the distribution
dist_sorted = sort(dist);

% Compute the 98th percentile
percentile_98 = dist_sorted(round(0.98 * length(dist)));

data = randi([14,2345],10000,1); % generate random integers between 14 and 2345
figure(2), clf
hist(data,100); % create a histogram with 100 bins
xlabel('Value'); % label the x-axis
ylabel('Frequency'); % label the y-axis

%% MATLAB for Brain and Cognitive Scientists (MIT Press 2017)
% mikexcohen@gmail.com
%
% This script contains solutions and tips for most exercises in the book.
% If you are struggling with an exercise not in this list, email me with an
% explanation of what you've already tried and what you need help with.

d = round(rand(10000,1)*(2345-14)+14);
figure(3), clf
histogram(d,50)

dd = prctile(d,98);

dsort = sort(d);
dsort(round(.98*length(d)));


title([ 'Range: ' num2str(min(d)) ' - ' num2str(max(d)) '; 98%: ' num2str(dd) ])


