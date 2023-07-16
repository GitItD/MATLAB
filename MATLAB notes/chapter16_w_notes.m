%% Matlab for Brain and Cognitive Scientists
% Matlab code for Chapter 16
% Mike X Cohen
% 
% This code accompanies the book 
%      "Matlab for Brain and Cognitive Scientists" (MIT Press, 2017).
% Using the code without following the book may lead to confusion, 
% incorrect data analyses, and misinterpretations of results. 
% Mike X Cohen assumes no responsibility for incorrect use of this code. 

%% generate bivariate correlated data
%{
This code generates a bivariate dataset with a sample size of 100, where the two variables are correlated with a 
Pearson correlation coefficient of 0.6. The first column of the data is filled with random values drawn from a 
standard normal distribution, while the second column is filled with the same values as the first column, 
transformed to have a correlation of 0.6 with the first column. The code then computes the covariance 
between the first variable in the dataset and a constant value of 100. Finally, the code creates a 
scatterplot of the two variables and calculates the observed Pearson correlation coefficient 
based on the generated data, which is displayed in the plot title.
%}
r = .6;
n = 100;
x = randn(n,2);
y = 100;
x(:,2) = x(:,1)*r + x(:,2)*sqrt(1-r^2);
%{
The line above updates the second column of the x matrix to create a 
correlation of r with the first column, while ensuring that the standard 
deviation of the second column remains the same. The updated values 
in the second column are calculated by taking the first column values (x(:,1)) 
and multiplying them by r, and then taking the existing values in the second column (x(:,2)) 
and multiplying them by the square root of (1-r^2). The resulting values in the 
second column have a correlation of r with the first column, as desired.
%}
cov1 = bsxfun(@minus, x, mean(x))'*(y - mean(y))/(n-1);     % Covariance formula

figure(1), clf
plot(x(:,1),x(:,2),'o')
obsR = corr(x);
title([ 'R=' num2str(obsR(2)) ])

%% multivariate covariance with real data
%{
This code loads in an EEG data and extracts a channels-by-time matrix from the first trial. 
It then calculates the covariance matrix of the data using two different methods: one where 
the mean of the data is not subtracted (the variable cov2), and another where the mean is 
subtracted (also stored in the variable cov2). The latter method is preferred because it calculates 
the covariance matrix of the deviation scores, which have a mean of zero.

The code then creates some visualizations of the covariance matrix and variances of the data. 
The first visualization displays the entire covariance matrix using an image plot. 
The second visualization shows the variances of the data along the diagonal of the 
covariance matrix using a topographic plot. The third visualization shows 
three seeded covariances (i.e., the covariances between specific electrode pairs) 
using topographic plots.
%}
load(['C:\Users\dun81\OneDrive\Documents\Dun II\ET Online\Library Online\Psychology\' ...
    'Computational Neuroscience\MATLAB for Brain and Cognitive Scientists - The MIT Press (2017)\MATLAB_BrainAndCognitiveScientists-main\sampleEEGdata.mat'])

% extract channels-by-time matrix from trial 1
data = squeeze(EEG.data(:,:,1));
%{
EEG.data is a 3D matrix that contains EEG data for multiple trials, 
multiple channels, and multiple time points. Each trial is a 3D slice 
along the third dimension of the matrix. squeeze is a MATLAB function 
that removes dimensions with only one element from a matrix, in this case, 
the third dimension of EEG.data. The resulting data matrix is a channels-by-time 
matrix from the first trial of the EEG data.
%}
% is the a valid covariance?
cov2 = data*data'/(EEG.pnts-1);
%{
This line of code calculates the sample covariance matrix for the multivariate 
data in the data variable. Specifically, it computes the covariance matrix as 
the product of data and its transpose data', divided by EEG.pnts-1.
EEG.pnts is one of the property of the data from EEG. 
%}

% this one's better (why?)
data = bsxfun(@minus, data, mean(data, 2));
cov2 = data*data'/(EEG.pnts-1);
%{
We want to minus the data by mean of time. Data of time is in the second dimension.
The bsxfun function is used to subtract the mean of the data along the second dimension 
from the data matrix. This is done to ensure that the mean of each channel is zero. 
We do this because when we calculate the covariance matrix, we want to measure 
the co-variation between channels and not their absolute values.

When we subtract the mean from a set of data, we are centering the data around zero. 
This means that the positive and negative values will be treated symmetrically, and the mean 
will no longer act as a measure of the magnitude of the data.

Normalizing the mean to zero has the effect of centering the data, which is important for 
computing covariance. The covariance between two variables X and Y is defined as the expected 
value of their product minus the product of their expected values:
%}

% show the covariance matrix
figure(2), clf
imagesc(cov2)
set(gca, 'clim', [-100 200])
%{
imagesc is a MATLAB function that displays a scaled image of a matrix. 
In this case, cov2 is a matrix containing the covariance values between 
all pairs of electrodes, which are estimated from the recorded EEG data. 
imagesc(cov2) displays a color-coded image of cov2, where each entry in the matrix is 
represented by a color, and the color scale is shown on the right-hand side of the image. 
The set(gca, 'clim', [-100 200]) command sets the color limits for the image to be between -100 and 200.
%}

% now show the variance ("co"variance along the diagonal)
%{
This code creates a figure with a topographic map of the variance values of the EEG data. 
The topoplotIndie function is used to plot the diagonal values of the cov2 matrix (i.e., the variances) 
on the scalp layout defined by the EEG.chanlocs structure. EEG.chanlocs is a field or property of a 
variable or structure EEG.

The set function is used to set the color limits of the plot to [0 200] and the colormap function is 
used to set the color scheme to 'hot'. The colorbar function is used to add a color bar to the 
plot and the title function is used to set the title of the figure to 'Variances'.
%}
figure(3), clf
topoplotIndie(diag(cov2), EEG.chanlocs);
set(gca,'clim',[0 200])
colormap hot
colorbar
title('Variances')

% and some seeded covariances. The numbers in cov2
% are particular electrodes.
%{
This code creates a figure with three subplots, each displaying a topographic 
plot of the covariance values between a specific EEG channel and all other channels.

In the first subplot, topoplotIndie is used to display the covariance values between 
the 47th row of cov2, using the corresponding electrode
 locations stored in EEG.chanlocs.
Similarly, the second and third subplots show the covariance values between 
the 23rd and 10th EEG channels, respectively, and all other channels.

The purpose of this code is to visualize the covariance values between each EEG channel 
and all other channels using topographic plots, which can help identify patterns and relationships between different channels.
%}
figure(4), clf
subplot(131)
topoplotIndie(cov2(47,:), EEG.chanlocs);

subplot(132)
topoplotIndie(cov2(23,:), EEG.chanlocs);

subplot(133)
topoplotIndie(cov2(10,:), EEG.chanlocs);

%% correlation in one variable
%{
This code generates two columns of random data, and uses them to 
demonstrate the difference between covariance and correlation. 
The variable r sets the correlation between the two columns. 
The code first calculates the covariance matrix using the cov function in Matlab. 
It then calculates the correlation using the formula cx / sqrt(vx1*vx2), 
where cx is the covariance, and vx1 and vx2 are the variances of the two columns. 
Finally, it calculates the correlation matrix using the corr function in Matlab. 
The two correlation values should be identical, but the covariance matrix will look 
different from the correlation matrix because the covariance is affected by
the scale of the data, while the correlation is not.
%}
r = .6;
n = 100;
x = randn(n,2);
x(:,2) = x(:,1)*r + x(:,2)*sqrt(1-r^2);
%{
This line of code creates a linear relationship between the first column of x 
and the second column of x, with correlation coefficient r. It generates a 
new value for the second column of x based on the values of the first 
column of x and some random noise. Specifically, it multiplies the 
first column of x by r, adds the second column of x multiplied by sqrt(1-r^2) 
(which is the square root of 1 minus the square of r) multiplied by some random 
noise, and stores the result in the second column of x. This produces a dataset 
in which the first and second columns have a linear relationship with correlation 
coefficient r.
%}

% stretch to make covariance different from correlation
x(:,1) = x(:,1)*100 + 2;
x(:,2) = x(:,2)*20 + 10;
%{
These two lines of code stretch the range of values in the first and second 
columns of x. The first column is multiplied by 100 and 2 is added to every 
element, while the second column is multiplied by 20 and 10 is added to every 
element. This is done to create a larger separation between the ranges of the 
two variables, making the difference between covariance and correlation more apparent.
%}
cx  = x(:,1)'*x(:,2) / (n-1);
vx1 = x(:,1)'*x(:,1) / (n-1);
vx2 = x(:,2)'*x(:,2) / (n-1);
%{
In the code, cx calculates the covariance between the first and second column of x, 
while vx1 and vx2 calculate the variances of the first and second column of x, respectively.

So the formulas for cx and vx are different because they measure different things. cx measures the joint 
variability of two variables, while vx measures the variability of each variable separately.
%}
cov1 = cov(x)

% these two variables should be identical!
cor1 = cx / sqrt(vx1*vx2)
cor2 = corr(x)
%{
cor1 and cor2 are two different ways to calculate the correlation between 
the two variables in matrix x. cor1 is calculated using the formula:

cor1 = cx / sqrt(vx1*vx2)

where cx is the covariance between the two variables, vx1 is the variance of the 
first variable, and vx2 is the variance of the second variable. This formula 
calculates the correlation coefficient using the definition of covariance and variance.

cor2 is calculated using the built-in MATLAB function corr, which calculates the 
correlation coefficient between columns of a matrix. corr(x) returns a matrix 
of correlation coefficients between all pairs of columns in x. In this case, since x 
has only two columns, corr(x) returns a 2x2 matrix with the correlation coefficient 
between the two columns in the upper-right and lower-left elements.

Both cor1 and cor2 should give the same result, as they are both calculating the 
correlation between the same two variables in x.
%}
%% correlation
%{
This code computes the correlation matrix of EEG data. Here are the steps:

1. The data matrix is extracted from the EEG object, with each row representing a channel and each column representing a time point.
2. The data is centered by subtracting the mean of each row from the corresponding row.
3. The data is normalized by dividing each row by its standard deviation.
4. The correlation matrix is computed as the covariance matrix of the normalized data.
5. The correlation matrix is visualized as an image using the imagesc function, with the color scale representing the correlation coefficients. 
The color range can be displayed using the get function applied to the clim (color limits) property of the current axis.
%}

% extract and prepare the data
data = squeeze(EEG.data(:,:,1));
data = bsxfun(@minus,data,mean(data,2));
data = bsxfun(@rdivide,data,std(data,[],2));


% compute correlation as covariance of normalized data
cor2 = data*data'/(EEG.pnts-1);

% image it
figure(5), clf
imagesc(cor2)
get(gca,'clim') % let's see the color range

%% or...
%{
The difference between the code above and below is how the covariance matrix is scaled. In the first code (above), 
the data is normalized by subtracting the mean and dividing by the standard deviation. 
Then the covariance matrix is computed from the normalized data, which gives the correlation matrix. 

In the second code (below), the covariance matrix is computed directly from the original data, without normalization. 
Then the variance matrix is computed from the diagonal of the covariance matrix, which is used to scale the 
covariance matrix. The scaled covariance matrix is displayed using imagesc.

The two methods differ in how they handle the scaling of the covariance matrix. 
The first method produces a correlation matrix, which has diagonal entries equal to 1 
and off-diagonal entries scaled by the standard deviations of the corresponding variables. 
The second method scales the covariance matrix by the product of the standard deviations of 
the corresponding variables, which makes the diagonal entries equal to the variances of the variables.
%}
% data again
data = squeeze(EEG.data(:,:,1));
data = bsxfun(@minus,data,mean(data,2));

% this is the unscaled covariance
cor2 = data*data'/(EEG.pnts-1);

% compute the variance matrix
stdMat = sqrt( diag(cor2)*diag(cor2)' );

% and then scale the (previously unscaled) covariance
cor2 = cor2 ./ stdMat;

imagesc(cor2)
get(gca,'clim') % how does this color range compare to that of the previous cell?

%% Anscobe's quartet to illustrate Pearson vs. Spearman
%{
The Anscombe's quartet is designed to demonstrate how the same basic statistical 
properties can be visually different, making it important to not only rely on numerical statistics.

In this specific case, each dataset has the same mean and variance, as well as the same 
correlation between the two variables. However, the data itself is generated differently, 
which causes it to have different distributions and therefore different shapes when plotted.

The first dataset is linear, the second dataset is non-linear, the third dataset has a single outlier, 
and the fourth dataset has one group of points with a different distribution. Therefore, when plotted, 
they will look different from each other, despite having similar statistical properties.

Anscombe's quartet is a set of four datasets that have the same mean, variance, correlation coefficient, 
and linear regression line, but they have different distributions and different relationships between the 
variables. This quartet is often used to illustrate the limitations of using only the Pearson correlation 
coefficient to measure the relationship between two variables.

The first dataset in the quartet has a linear relationship between the variables, so the Pearson correlation 
coefficient accurately represents the strength of the relationship. However, the other three datasets have 
non-linear relationships between the variables, and in these cases, the Pearson correlation coefficient is 
not a good measure of the strength of the relationship.

The second dataset has a strong non-linear relationship, and the Pearson correlation coefficient suggests 
that the relationship is weak, which is misleading. The third dataset has an outlier that has a large effect on 
the correlation coefficient, but the Spearman rank correlation coefficient is not affected by outliers, and it 
accurately represents the strength of the relationship. The fourth dataset has a relationship that is entirely 
dependent on a single point, which has a large effect on the Pearson correlation coefficient, but again, the 
Spearman rank correlation coefficient is not affected by this issue.

So, by looking at the four datasets in Anscombe's quartet, one can see that the Spearman rank correlation 
coefficient is a more robust measure of the strength of the relationship between two variables when the 
relationship is non-linear or there are outliers.
%}
anscombe = [
  % series 1    series 2    series 3     series 4
    10  8.04    10  9.14    10  7.46     8  6.58;
     8  6.95     8  8.14     8  6.77     8  5.76;
    13  7.58    13  8.76    13 12.74     8  7.71;
     9  8.81     9  8.77     9  7.11     8  8.84;
    11  8.33    11  9.26    11  7.81     8  8.47;
    14  9.96    14  8.10    14  8.84     8  7.04;
     6  7.24     6  6.13     6  6.08     8  5.25;
     4  4.26     4  3.10     4  5.39     8  5.56;
    12 10.84    12  9.13    12  8.15     8  7.91;
     7  4.82     7  7.26     7  6.42     8  6.89;
     5  5.68     5  4.74     5  5.73    19 12.50;
    ];


%% end
