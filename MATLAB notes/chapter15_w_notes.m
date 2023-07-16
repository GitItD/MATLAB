%% Matlab for Brain and Cognitive Scientists
% Matlab code for Chapter 15
% Mike X Cohen
% 
% This code accompanies the book 
%      "Matlab for Brain and Cognitive Scientists" (MIT Press, 2017).
% Using the code without following the book may lead to confusion, 
% incorrect data analyses, and misinterpretations of results. 
% Mike X Cohen assumes no responsibility for incorrect use of this code. 

%% theoretical and empirical distributions
%{
This code creates two subplots in a single figure window (figure(1)). 
In the first subplot, it plots the normal distribution and marks the observed value 
with a vertical line in magenta color. In the second subplot, it plots the histogram 
of the permuted mean height differences and marks the observed value with a vertical line in magenta color.

The permutation test is used to test the null hypothesis that there is no difference between the two groups 
in the heights. The mean height difference between the two groups is calculated and stored in the 
"observed_diff" variable. The permutation test is performed by randomly shuffling the group assignment 
1000 times and recalculating the mean height difference for each shuffle. The resulting mean height 
differences are stored in the "permdiffs" array.

The purpose of this code is to show the difference between the "theoretical" distribution and 
the "empirical" distribution by plotting both in separate subplots.
%}
figure(1), clf

obsval = 1.8;
% Load sample data
heights = [1.75, 1.65, 1.8, 1.72, 1.68, 1.9, 1.75, 1.82, 1.73, 1.85];
group_assignment = [1, 1, 2, 2, 2, 1, 1, 2, 2, 2];

% Calculate observed mean height difference
%{
The code calculates the difference in mean heights between two groups. 
The heights data is divided into two groups based on the values in the "group_assignment" variable. 
The mean of the heights for group 1 is calculated and then subtracted 
from the mean of the heights for group 2, giving the observed_diff value.

In this code, the groups 1 and 2 appear to be assigned in an intentional manner, as the heights of the 
individuals are separated based on their group assignment and the mean height is calculated for 
each group. However, it is not specified in the code how the group assignment was 
determined, so it is possible that the assignment was random or based on some other criterion.
%}
observed_diff = mean(heights(group_assignment == 1)) - mean(heights(group_assignment == 2));

% Set the number of permutations
%{
The number of permutations is the number of times the code will repeat the random 
shuffling of the group assignments. In this case, the number of permutations is set to 
1000, so the random shuffling will be done 1000 times.
%}
num_permutations = 1000;

% Initialize an array to store the permuted mean height differences
%{
This line initializes an array called "permdiffs" with 1000 rows and 1 column, and 
sets all the values to zero. The purpose of this array is to store the mean 
height differences for each permutation.
%}
permdiffs = zeros(num_permutations, 1);

% Perform permutation test
%{
This code performs a permutation test by shuffling the group assignments num_permutations 
times and calculating the mean height difference for each shuffled group assignment. 
The mean height difference is calculated by subtracting the mean of the heights in 
group 2 from the mean of the heights in group 1. The result of each 
permutation is stored in the permdiffs array.

Yes, the for loop that performs the permutation test is the most important part of 
this section as it calculates the mean height differences for the 1000 permuted samples. 
This generates the empirical distribution of the mean height differences, which can then be 
compared to the observed mean height difference to determine the significance of the 
difference between the two groups

randperm is a function in MATLAB that returns a random permutation of the integers from 1 to n, 
where n is a positive integer specified as an input argument. 

The line shuffled_assignment = group_assignment(randperm(length(group_assignment))) creates a 
new variable shuffled_assignment that is a randomly permuted version of the original group_assignment.

This line of code calculates the mean height difference between the two groups after shuffling the assignment 
of the individuals to the groups. The mean height of each group after shuffling is found using mean(heights(shuffled_assignment == 1)) 
and mean(heights(shuffled_assignment == 2)). The difference of these two values gives the mean height difference 
between the two groups after shuffling, which is then stored in permdiffs(i). 
The loop runs num_permutations times, so permdiffs will store num_permutations mean 
height differences after shuffling the assignment of individuals to the groups.
%}
for i = 1:num_permutations
    shuffled_assignment = group_assignment(randperm(length(group_assignment)));
    permdiffs(i) = mean(heights(shuffled_assignment == 1)) - mean(heights(shuffled_assignment == 2));
end

% plot "theoretical" distribution
subplot(121)
x=-3:.001:3;
plot(x,exp(-x.^2))
hold on
plot([obsval, obsval],get(gca,'ylim'),'m')
set(gca,'xlim',[-4 4])
xlabel('Statistic values')
ylabel('Probability')


% plot empirical distribution
% NOTE: You must first create Figure 15.3 for the next few lines to work.
subplot(122)
histogram(permdiffs, 50) 
hold on
plot([obsval obsval],get(gca,'ylim'),'m')
set(gca,'xlim',[-4 4])
xlabel('Statistic values')
ylabel('Count')

%% shifted and slightly skewed distribution
%{
This code generates two histograms, one of a Gaussian distribution (generated by randn(10000,1)) and 
one of a skewed distribution. The Gaussian distribution is first skewed by transforming the values 
that are greater than 0 with the logarithm log(1+r(r>0)). 
This transformation makes the distribution slightly more positive-skewed. 
After skewing the distribution, the values are shifted to make them all positive by 
subtracting the minimum value from all values and negating the result. 
The two histograms are then plotted and displayed with a legend indicating which is which.
%}

% create Gaussian distribution
r = randn(10000,1);

% get histogram values for later comparison
%{
The line [y1,x1] = hist(r,100); computes the histogram values of the Gaussian distribution r with 100 bins. 
The output y1 gives the number of elements in each bin, and x1 gives the bin edges. 
This histogram can be thought of as a rough estimate of the probability density function (pdf) of the distribution.
%}
[y1,x1] = hist(r,100);


% skew the distribution
%{
This code is transforming the Gaussian distribution to create a positive skewed distribution. 
The logarithmic transformation is applied to all positive values of the Gaussian distribution. 
By adding 1 to the argument of the log function, the function is guaranteed to be defined 
for all positive values. The resulting distribution is then shifted to the right to make it all positive 
by subtracting the minimum value.
%}
r(r>0) = log(1+r(r>0));
% and make it positive
r = -r-min(r);

% and get its histogram values
[y2,x2] = hist(r,100);



figure(2), clf
plot(x1,y1), hold on
plot(x2,y2,'r')
legend({'Gaussian';'positive skewed'})


%% simulated firing rate distributions
%{
Similar as the codes right above,
%}
% number of trials
N = 100;

% male pictures
%{
The line fr_males = 26-r*10 calculates the firing rate (fr_males) for male pictures. 
It starts by creating a normally distributed random variable (r = randn(N,1)), 
then transforms the positive values of the distribution by applying the logarithmic function (r(r>0) = log(1+r(r>0))). 
Finally, the firing rate for male pictures is calculated by subtracting the transformed random variable (r) from 26 
and multiplying the result by 10 (fr_males = 26-r*10).

The value of 26 here is an arbitrary value used to simulate firing rates for the "male" group in this code. 
It represents the mean firing rate for the "male" group, but you could use any other value. 
The standard deviation of the firing rates for this group is controlled by the value of 10, which is used to scale the 
random numbers generated from the normal distribution. A smaller value of 10 would result in a smaller standard 
deviation and tighter clustering of the firing rates, while a larger value of 10 would result in a 
larger standard deviation and more spread out firing rates.
%}
r = randn(N,1);
r(r>0) = log(1+r(r>0));
fr_males = 26-r*10;

% get histogram values for later comparison
[y1,x1] = hist(fr_males,20);


% female pictures
r = randn(N,1);
r(r>0) = log(1+r(r>0));
fr_females = 30-r*10;

% get histogram values for later comparison
[y2,x2] = hist(fr_females,20);


figure(3), clf
plot(x1,y1), hold on
plot(x2,y2,'r')
legend({'Males';'Females'})

%% mix trials together
%{
This code concatenates two arrays, fr_males and fr_females, into a single array allfr. 
Then, it mixes the trials up by generating a random permutation of the indices of allfr 
using the randperm function and reordering the elements of allfr based on the permuted indices. 
Finally, the conditions corresponding to the trials are relabeled as 0 or 1 based on their order in the permuted allfr array.
%}

% concatenate trials
allfr = cat(1,fr_males,fr_females);

% the following line does the same thing as the previous
allfr = [fr_males, fr_females];

% and mix them up (but in what order??)
allfr = allfr(randperm(N*2));
% a better approach
allfr = [fr_males, fr_females];
% repeat original and new order variables for accurate tracking
%{
The deal function creates multiple output variables, but each one is assigned the same input, which is randperm(N*2). 
In this case, the output variables are seperated into conds and neworder.
%}
[conds, neworder] = deal( randperm(N*2) );

% re-sort data
%{
This line sorts the allfr variable based on the order defined in neworder. 
The allfr vector is reordered so that the elements are in the order specified by the neworder vector.
%}
allfr = allfr(neworder);

% and relable the conditions
%{
This line of code takes the neworder vector, which is a permutation of the numbers 1 to N*2, 
and replaces all the numbers in it that are less than N+1 with the value 1. 
All other numbers greater than N+1 are then replaced with the value 0. 
So this line of code is essentially re-labeling the conditions based on the new permutation.
%}
conds(neworder<N+1) = 1;
conds(conds>1) = 0;

%% generate one null hypothesis scenario
%{
The outputs of these two sets of mean calculations will show if the means are different or not. 
If the means are different, it suggests that the original conditions are not randomly assigned 
and there is a significant difference between the firing rates for males and females. 
However, if the means of the fake conditions are not different, it suggests that the 
null hypothesis that the firing rates for males and females are not significantly different can not be rejected.
%}

% random permutation
fakeconds = randperm(N*2);

% shuffled condition labels
fakeconds(fakeconds<N+1) = 1;
fakeconds(fakeconds>1) = 0;


% these two means should be different.
[mean(allfr(conds==1)) mean(allfr(conds==0))]

% should these two be different?
[mean(allfr(fakeconds==1))  mean(allfr(fakeconds==0)) ]

%% and now a distribution of null hypothesis values
%{
This code generates 1000 random permutations of the trial labels and calculates the mean firing rate 
difference between the two conditions (males and females) for each permutation. 
The distribution of these mean differences under the null hypothesis (that the firing 
rate difference between the two conditions is random) is plotted in a histogram using the hist function. 
The observed mean firing rate difference between the two conditions is then plotted as a vertical 
line on top of the histogram using the plot function. The x-axis represents the value of the mean firing rate 
difference and the y-axis represents the count of how many times that value appears in the distribution.
%}
nPerms = 1000;
permdiffs = zeros(nPerms,1);

for permi=1:nPerms
    fconds = randperm(N*2);
    fconds(fconds<N+1) = 1;
    fconds(fconds>1) = 0;
    permdiffs(permi) = mean(allfr(fconds==0))-mean(allfr(fconds==1));
end


% plot the distribution of H0 values
figure(4), clf
hist(permdiffs,50)
hold on

% and plot the observed value on top
obsval = mean(allfr(conds==0))-mean(allfr(conds==1));
plot([obsval obsval],get(gca,'ylim'),'m','linew',10)
xlabel('Value'), ylabel('Count')

%% two methods of evaluating statistical significance
%{
The first method, the Z-value method, uses the Z-score to calculate the statistical significance of the observed 
difference in means. The Z-score is the difference between the observed value and the mean of the distribution 
of differences under the null hypothesis, divided by the standard deviation of that distribution. 
This Z-score is then transformed into a p-value, which is the probability of observing a difference in 
means as extreme as the observed difference under the null hypothesis.

The second method, the p-value count method, simply counts the number of times the difference in means 
under the null hypothesis is greater than the observed difference in means, and divides by the total number 
of permutations. This gives an estimate of the p-value, the probability of observing a difference in means as 
extreme or greater than the observed difference under the null hypothesis.
%}

% Z-value
zVal = ( obsval-mean(permdiffs) ) / std(permdiffs);
p = 1-normcdf(abs(zVal));

% p-value count
pCount = sum(permdiffs>obsval)/nPerms;

%% now with real data
%{
This code is performing a two-sample t-test on the power difference of two conditions 
(presumably some kind of experiment with two conditions). The first 100 timepoints are 
assigned to one condition and the next 100 timepoints are assigned to the second condition. 
The difference in power between the two conditions is computed and stored in realdif.

The permutation testing is then performed by randomly reassigning the 100 timepoints 
to the two conditions 1000 times, computing the difference in power for each permutation, 
and storing the results in permdif.

Next, a z-score is calculated for each time-frequency point using realdif and the mean 
and standard deviation of permdif. This z-score map is stored in zmap.

Finally, the code plots four subplots:
1) The power difference between the two conditions.
2) The z-score map.
3) A thresholded z-score map, where all z-scores below a certain threshold (determined by the p-value of .05) are set to 0.
4) The power difference between the two conditions, with significant regions (determined by the thresholded z-score map) outlined with black contours.
%}
load tfdata

% p-value threshold
pval = .05;

% the empirically observed time-frequency power difference
%{
squeeze is a MATLAB function that removes singleton dimensions from an array. 
For example, if A is a matrix (a 2-dimensional array), squeeze(A) will return the same matrix. 
However, if A is a 3-dimensional array with a singleton dimension (i.e., one of the dimensions has size 1), 
squeeze(A) will return a 2-dimensional array with that singleton dimension removed.

In the code snippet you provided, allpow is a 3-dimensional array with dimensions channel x frequency x time. 
The line realdif = squeeze( mean(allpow(:,:,101:200),3) - mean(allpow(:,:,1:100),3) ); calculates the difference 
between the mean power in the first 100 time points and the mean power in the next 100 time points 
for each channel and frequency, resulting in a 2-dimensional array with dimensions channel x frequency. 
The squeeze function is used to remove any singleton dimensions that may have been created 
during the calculation of the mean power.
%}
realdif = squeeze( mean(allpow(:,:,101:200),3) - mean(allpow(:,:,1:100),3) );

% setup permutation testing
nPerms  = 1000; % number of iterations
permdif = zeros([ nPerms size(realdif) ]);

% loop over permutations
for permi=1:nPerms
    
    % generate a fake trial order
    fakeord = randperm(size(allpow,3));
    
    % compare this to how the variable 'realdif' was created
    fakedif = squeeze( mean(allpow(:,:,fakeord(1:100)),3) - mean(allpow(:,:,fakeord(101:200)),3) );
    
    % does this need to be a separate line from the previous?
    permdif(permi,:,:) = fakedif;
end

% compute z difference score
%{
The "mean(permdif, 1)" function calculates the mean along dimension 1 of the "permdif" matrix. 
The result is a row vector containing the mean of each column of the "permdif" matrix.
In this case, the first dimension refers to the rows of the permdif matrix. The syntax mean(permdif, 1) 
calculates the mean of each row along the first dimension of the matrix. 

This code computes the standard deviation along the first dimension of the "permdif" matrix. 
The result will be a 2D matrix with the same size as the second and third dimensions of the "permdif" matrix, 
with each element in the matrix representing the standard deviation of the values along the first dimension 
for that particular position. 

The second argument, "[]" tells the "std" function to compute the standard deviation over 
all the elements along that dimension, so the result will be a single number for each position in the 2D matrix. 
The "squeeze" function will remove any singleton dimensions from the result.
%}
permmean = squeeze(mean(permdif,1));
permstd  = squeeze(std(permdif,[],1));
zmap = ( realdif-permmean ) ./ permstd;


% plotting...
%{
The code colorval = max(abs(realdif(:)))*.8 calculates the maximum absolute value of 
realdif and multiplies it by 0.8. The purpose of this line of code is likely to determine a 
color scale for visualizing realdif later in the code. The value of colorval would be used 
as a threshold for the color scale, such that values greater than colorval will be shown 
as the maximum color, and values less than -colorval will be shown as the minimum color.
The operator (:) is used to convert a multi-dimensional matrix into a one-dimensional vector.
%}
figure(5), clf
subplot(221)
contourf(timevec, frex, realdif, 40, 'linecolor', 'none')
colorval = max(abs(realdif(:)))*.8;
set(gca, 'clim',[-colorval, colorval])
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
title('Power condition difference')


subplot(222)
contourf(timevec, frex, zmap, 40, 'linecolor', 'none')
colorval = max(abs(zmap(:)))*.8;
set(gca,'clim', [-colorval, colorval])
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
title('Power condition Z-map')


subplot(223)

% Threshold the map by turning all subtreshold pixels to 0.
zthresh = zmap;
zthresh(abs(zthresh)<norminv(1-pval)) = 0;

contourf(timevec,frex,zthresh,40,'linecolor','none')
colorval = max(abs(zthresh(:)))*.8;
set(gca,'clim',[-colorval colorval])
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
title('Thresholded differences')

subplot(224)
contourf(timevec,frex,realdif,40,'linecolor','none')
colorval = max(abs(realdif(:)))*.8;
set(gca,'clim',[-colorval colorval])
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
title('Significant regions outlined')

% Draw significance contours on top of the colorful map.
hold on
contour(timevec,frex,logical(zthresh),1, 'k')

%% extreme-value-based correction
%{
This code performs an extreme-value-based correction for multiple comparisons. 
It first creates a random permutation of the trials (fakeord) and calculates the difference 
in power between the first 100 trials and the next 100 trials for each permutation (fakedif). 
Then it finds the most extreme positive and negative values for each permutation and 
stores these in the exvals matrix. The code then calculates the lower and upper thresholds 
using the prctile function, based on the desired p-value (here set to 0.05). 
Finally, the code plots a histogram of the distribution of the extreme values 
and superimposes a histogram of the real values, along with the calculated 
lower and upper thresholds.
%}

% p-value threshold
pval = .05;

% initialize extreme-value distribution
exvals = zeros(nPerms, 2);

for permi=1:nPerms
    
    % again with the shuffled trial order
    %{
    fakeord is a 1D array of indices, so if allpow is a 3D array, then allpow(:,:, fakeord) 
will shuffle the third dimension of allpow according to the indices specified in fakeord.
    %}
    fakeord = randperm(size(allpow, 3));    
    fakedif = squeeze(mean(allpow(:,:, fakeord(1:100)), 3) - mean(allpow(:,:, fakeord(101:200)), 3) );

    % from this map, take the most extreme positive and negative values
    exvals(permi, 1) = min(fakedif(:));
    exvals(permi, 2) = max(fakedif(:));
end

% find the thresholds
%{
The variable lowerThresh is assigned the 100*pval percentile value of the first 
column (column 1) of the exvals matrix, while the variable upperThresh is assigned 
the 100-100*pval percentile value of the second column (column 2) of the exvals matrix. 
The prctile function is used to calculate percentiles, which are values that divide a set of 
observations into 100 equal parts, with each part representing 1 percentile.
%}
lowerThresh = prctile(exvals(:,1), 100*pval);
upperThresh = prctile(exvals(:,2), 100 - 100*pval);


% let's have a look at the distribution of extreme values
figure(6), clf
histogram(exvals(:), 200, 'DisplayStyle','stairs')
xlabel('Value'), ylabel('Count')

% plot the distribution of real values on top
hold on
h=histogram(realdif(:), 200);

% and draw the thresholds
plot(ones(1, 2)*lowerThresh, get(gca,'ylim'), 'k--')
plot(ones(1, 2)*upperThresh, get(gca,'ylim'), 'k--')

% finally, change some of the histogram properties
set(h,'edgecolor','none','facecolor','r')

%%

% threshold the real map
threshmap = realdif;
threshmap(threshmap>lowerThresh & threshmap<upperThresh) = 0;



figure(7), clf
subplot(221)
contourf(timevec,frex,realdif,40,'linecolor','none')
colorval = max(abs(realdif(:)))*.8;
set(gca,'clim',[-colorval colorval])

subplot(222)
contourf(timevec,frex,threshmap,40,'linecolor','none')
colorval = max(abs(threshmap(:)))*.8;
clim_vals = [-colorval colorval];
if all(isnumeric(clim_vals)) && length(clim_vals) == 2 && clim_vals(2) > clim_vals(1)
    set(gca, 'CLim', clim_vals)
else
    disp('Error: Invalid value for clim_vals')
end


subplot(223)
contourf(timevec,frex,realdif,40,'linecolor','none')
colorval = max(abs(realdif(:)))*.8;
set(gca,'clim',[-colorval colorval])

hold on
contour(timevec,frex,logical(threshmap),1,'k')

%% end
