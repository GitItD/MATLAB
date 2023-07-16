%% Matlab for Brain and Cognitive Scientists
% Matlab code for Chapter 14
% Mike X Cohen
% 
% This code accompanies the book 
%      "Matlab for Brain and Cognitive Scientists" (MIT Press, 2017).
% Using the code without following the book may lead to confusion, 
% incorrect data analyses, and misinterpretations of results. 
% Mike X Cohen assumes no responsibility for incorrect use of this code. 

%{
	1. Import and pre-process data: The raw data needs to be imported and pre-processed to ensure that it meets the requirements of signal detection theory.
	2. Define stimulus levels: The different stimulus levels must be defined and assigned to the data.
	3. Define response categories: The response categories must be defined, such as hit, miss, false alarm, and correct rejection.
	4. Compute hit rate and false alarm rate: These rates can be computed using MATLAB functions like find, bsxfun, and histc.
	5. Plot Receiver Operating Characteristic (ROC) curve: An ROC curve can be plotted using MATLAB functions such as plot and nbins.
	6. Compute conditional accuracy function: The conditional accuracy function can be computed using MATLAB functions such as linspace and binidx.
	7. Analyze results: The results of the analysis can be analyzed using MATLAB functions like tiedrank and caf.
Visualize results: The results can be visualized using MATLAB functions like dp2plot and ind2sub.
%}

%% Computing Discrimination d'

% step 1: convert hits and false alarms from counts to proportions;
hitP = 22/30;
faP  =  3/30;

% step 2: convert those proportions to standard z-values
hitZ = norminv(hitP);
faZ  = norminv(faP);

% step 3: take the difference of those z-values
dPrime = hitZ-faZ;

%% aside on norminv

x = 0:.001:1;

p1 = norminv(x); % with stats toolbox
p2 = -sqrt(2).*erfcinv(2*x); % without stats toolbox, also works in Octave
%{
erfcinv is the inverse complementary error function, which maps a probability 
in the interval [0, 2] to a standard normal deviate. It is used in mathematical and 
statistical calculations to find the inverse of the cumulative distribution function 
for the standard normal distribution. It's worth noting that the 
complementary error function only covers the upper half of the standard Gaussian curve.
%}

figure(1), clf
plot(p1,x,'b','linew',8), hold on
plot(p2,x,'ro')
xlabel('Z value'), ylabel('proportion')

%% 2D d' space
%{
The code below is used to plot a receiver operating characteristic (ROC) or isosensitivity curve in signal detection theory. 
The ROC plot is a graphical representation of the relationship between the hit rate and false alarm rate for a particular signal detection task. 
The dŒ value is a measure of pure discriminability, uncoupled from the overall rate of responding, which makes it robust to response biases.

To create the ROC plot, a 2D matrix of dŒ values is created with hit rate on the y-axis and false alarm rate on the x-axis. 
This is done by having the MATLAB function bsxfun expand a row vector and a column vector into a matrix and then
 perform point-wise subtraction between those matrices. This allows for the calculation of dŒ for many values 
simultaneously. The resulting ROC plot displays the ranges of hit and false alarm rates that produce the same sensitivity (dŒ).
%}

% convert counts to probabilities
x  = .01:.01:.99;

dp = bsxfun(@minus,norminv(x)',norminv(x));
%{
The code above creates a matrix of d' values, with hit rate on the y-axis and false alarm rate on the x-axis. 
The bsxfun function is used to point-wise subtract two matrices, which are created by expanding a row vector and a column vector.

The row vector norminv(x)' represents the normal inverse cumulative distribution function applied to the false alarm rate, 
while the column vector norminv(x) represents the normal inverse cumulative distribution function applied to the hit rate. 
The difference between the two is the d' value, which measures the difference between the means of two normal distributions 
representing the signal and noise in a signal detection experiment.
%}

% show the 2D d' space
figure(2), clf
contourf(x,x,dp,80,'linecolor','none')
xlabel('False alarm rate')
ylabel('Hit rate')
hold on % we'll plot lines on top
axis square
%{
The code above creates a contour plot of the discriminability index dŒ as a function of false 
alarm rate and hit rate. The contourf function creates a filled contour plot with 80 contour lines, 
and the linecolor option is set to 'none' so that the contour lines are not displayed. 
The x-axis is labeled as "False alarm rate" and the y-axis is labeled as "Hit rate". 
The hold on command is used so that additional lines can be plotted on top of the contour plot. 
The axis square command sets the aspect ratio of the plot to be equal, so that the x and y scales are the same.

contourf(x,x,dp,80,'linecolor','none') creates a filled contour plot of the 2D d' space. The parameters are as follows:
x: a row vector of false alarm rates, ranging from .01 to .99 in increments of .01
x: a row vector of hit rates, ranging from .01 to .99 in increments of .01
dp: a 2D matrix of d' values, computed as dp = bsxfun(@minus,norminv(x)',norminv(x))
80: the number of contour levels to plot
'linecolor','none': specifies that the lines separating the contour levels should not be plotted, only the filled contours.
The contourf function plots the 2D d' space by showing the values of d' as different levels of shading or color. 
The x-axis represents the false alarm rate and the y-axis represents the hit rate.
%}

% colors for isosensitivity curves
colorz = 'rbmk';
%{
The value "colorz" is a character array that defines a set of colors for plotting the isosensitivity curves. 
The letters 'r', 'b', 'm', and 'k' stand for red, blue, magenta, and black, respectively, which are standard color codes in MATLAB.
%}

% the d' values
dp2plot = [1 1.5 2 2.5];
tol = .01;
%{
The d' values used in the code are arbitrary values specified by the programmer: dp2plot = [1 1.5 2 2.5]. 
They are used to plot the isosensitivity curves for different d' values.
%}

for dpi=1:length(dp2plot)
    
    % find points
    idx = find(dp>dp2plot(dpi)-tol & dp<dp2plot(dpi)+tol);
    
    % and plot isosensitivity curves
    [yi,xi] = ind2sub(size(dp), idx);
    plot(x(xi),x(yi),[ colorz(dpi) 'o-' ],'linew',4,'markersize',9)
end
%{
This is MATLAB code to plot isosensitivity curves. 

The dp2plot variable specifies the dŒ values that you want to plot. 
The code starts a loop that will go through each dŒ value. 
For each dŒ value, it first finds the indices of the points in the dŒ
matrix (dp) that have dŒ values within tolerance range in the + or - of
 tol.

The ind2sub function converts these indexes from linear indexes to 2D indexes, which are the row and column indexes of the points. 
Finally, the code plots these points and connects them with a line, with the color and style specified by the colorz and '-o' options, respectively. 
The 'linew' and 'markersize' options control the width of the lines and the size of the markers, respectively.
%}

%% response bias
%{
The code below calculates dŒ and response bias from hit rate (hitP) and false alarm rate (faP).

First, it calculates the z-scores for hit rate and false alarm rate using the norminv function. 

Next, the code calculates dŒ as the difference between the z-scores of the hit rate and false alarm rate.

Finally, the code calculates the response bias, which is the average of the 
negative of the sum of the z-scores of hit rate and false alarm rate.
%}
% step 1
hitP = 22/30;
faP  =  3/30;

% step 2
hitZ = norminv(hitP);
faZ  = norminv(faP);

% step 3
dPrime = hitZ-faZ;
respBias = -(hitZ+faZ)/2;

%% 2D response bias space
%{
This code is generates a 2D plot of response bias in a signal detection task.
The script first calculates the response bias metric based on the hit rate 
(probability of detecting a signal when it is present) and false alarm rate 
(probability of detecting a signal when it is absent).

The script then creates a contour plot of the response bias using the hit 
and false alarm rates as the x and y axes, respectively. 
The plot is annotated with various "isosensitivity curves" which are sets of points 
with a constant response bias. These isosensitivity curves are plotted on top of the contour 
plot and are visualized as a series of points connected by a line.

Overall, this script provides a visual representation 
of the relationship between hit and false alarm rates 
and the resulting response bias in a signal detection task.
%}

% convert counts to probabilities
rb = -bsxfun(@plus,norminv(x)',norminv(x))/2;
%{
The code creates a 2D bias space. rb is a matrix of response bias values, 
and the calculation is performed using the bsxfun function in Matlab.

The bsxfun function performs binary operations element-wise on arrays of different sizes. 
In this case, it is subtracting the transpose of norminv(x) from norminv(x) and dividing the result by 2.

norminv is the inverse cumulative distribution function (CDF) of the standard normal distribution, 
which is used to convert probabilities to z-scores. The calculation of response bias involves subtracting 
the z-scores for the false alarm rate and hit rate, so the calculation is done using the norminv function. 
The division by 2 is a standard conversion for the calculation of response bias in signal detection theory.
%}

figure(3), clf
contourf(x,x,rb,80,'linecolor','none')
xlabel('False alarm rate')
ylabel('Hit rate')
hold on
axis square
%{
This code creates a plot of the "Response Bias Space". 
It plots the response bias values as a function of the false alarm rate and hit rate. 

The contourf function creates a filled contour plot of the response bias values, with 80 levels, 
and with no outlines around the contour levels. 

The xlabel and ylabel functions add labels to the x and y axis of the plot, respectively. 
The hold on function holds the plot so that subsequent plotting commands can be added to it. 
The axis square function sets the aspect ratio of the plot to be equal in both x and y dimensions.
%}
colorz = 'rbmk';

rb2plot = [.3 .5 .9 1.5];
tol = .01;

for dpi=1:length(rb2plot)
    
    % find points
    idx = find(rb>rb2plot(dpi)-tol & rb<rb2plot(dpi)+tol);
    
    % and plot isosensitivity curves
    [yi,xi] = ind2sub(size(rb),idx);
    plot(x(xi),x(yi),[ colorz(dpi) 'o-' ],'linew',4,'markersize',9)
end
%{
The above loops is the same as 2D d' space.
%}

%% discretization
%{
Discretization is the process of transforming a continuous variable into a categorical (discrete) variable. 
It involves dividing a continuous range into a number of "bins" and then representing the values 
within a bin by a single, representative value, such as the mean. 
The goal is to reduce the number of values of the continuous 
variable and make the data easier to analyze and interpret. 
This process is used in various fields such as data science, machine learning, and statistics.

Discretization is the process of dividing a continuous variable into discrete, equal-width or unequal-width categories, or bins. 
The steps of discretization are:

Determine the number of bins to use.
Sort the data in ascending order if using equal-width bins.
Divide the range of the data into equally spaced bins or into bins of unequal width based on the desired bin sizes.
Assign each data point to its nearest bin.
Compute the mean of each bin for each of the relevant variables.
Plot the discretized data and interpret the results.

To discretize a vector based on values in the data, we have two options:
the easy way and the better way. The easy way is to sort the data in ascending
order and then use the function linspace to define discretization bin
edges. In this example, data from 100 trials and two columns will be discretized.
%}
% define parameters
ntrials = 100;
nbins   = 7;

% create random data
d = [500+100*randn(ntrials,1) rand(ntrials,1)>.3];
d = sortrows(d,1);

% discretization based on n bins
binidx = ceil(linspace(0,nbins-1,length(d)));
discdata = zeros(2,nbins);

for i=1:nbins
    discdata(1,i) = mean(d(binidx==i,1));
    discdata(2,i) = mean(d(binidx==i,2));
end

% and plot.
figure(4), clf
plot(discdata(1,:),discdata(2,:),'o-','markerfacecolor','r','markersize',10)

%% discretization using tiedrank
%{
This code is performing a better way of discretizing data into N roughly equally sized groups.
 
The data is generated as random values with 100 trials and two columns. 

The discretization process is performed in two steps. 
The first step is to calculate the tiedrank of the first 
column of the data and normalize it by the number of trials (ntrials). 
The second step is to multiply the tiedrank by the number of bins (nbins) 
and round up to the nearest integer using the function ceil. 
The result of this process is stored in the variable drank.

This code is considered better because it allows the original 
data to remain in its original order, preserving trial order 
which might be important for other analyses. 

The tiedrank function is used to increase the precision of the discretization. 
The tiedrank function assigns ranks to the data based on their values, 
such that data values that are tied receive the same rank, 
and the mean rank is assigned to the tied data points. 

By using the tiedrank function, the resulting discretized data will have 
more precision compared to the previous method of using the linspace function.
%}

% more random data
d = [500+100*randn(ntrials,1) rand(ntrials,1)>.3];

% procedure in separate steps
temp  = tiedrank(d(:,1))/ntrials; % tiedrank is not in default Octave, but you can find it online
temp  = nbins*temp;
drank = ceil( temp );
%{
tiedrank: Rank adjusted for ties
Examples:
Counting from smallest to largest, the two 20 values are 2nd and 3rd, so they both get rank 2.5 (average of 2 and 3):

tiedrank([10 20 30 40 20])
ans =
    1.0000    2.5000    4.0000    5.0000    2.5000 
%}

% can also be done in one line
drank = ceil( nbins*tiedrank(d(:,1))/ntrials );

%% CAFs (Confidence Accuracy Function) with real data
%{

The Conditional Accuracy Function (CAF) is a relationship between performance accuracy and response time.

To calculate a CAF:
	1. Discretize the data into N roughly equally sized groups. This can be done by calculating the rank of 
each data point and dividing it by the total number of data points, then multiplying by N and rounding up to the 
nearest integer to determine the group that each data point belongs to.
	2. For each group, calculate the average performance of the data points in that group.
	3. Plot the average performance for each group on the Y-axis and the corresponding group midpoint on the X-axis.

The result is a plot of the CAF, showing the relationship between performance and the underlying variable being discretized.
%}

load behavioralDataRK.mat

% number of discretizations
nbins = 12;

% and now discretize
drank = ceil( nbins*tiedrank(beh(:,2))/length(beh) );
%{
The line of code calculates the discretized rank of the data in column 2 of the matrix "beh". 
The "tiedrank" function returns the fractional rank of each value in the data, and then this is 
multiplied by the number of bins (nbins) to stretch the range of values. 
The "ceil" function rounds up to the nearest integer, so the resulting drank 
variable is an integer indicating which bin each value in column 2 of "beh" belongs to.

In the code, "beh" refers to an input variable/matrix that has been loaded from a .mat file.
%}

% initialize and compute
caf = zeros(nbins,2);
for i=1:nbins
    caf(i,1) = mean(beh(drank==i,2));
    caf(i,2) = mean(beh(drank==i,1));
end
%{
This code calculates the Conditional Accuracy Function (CAF). 
It starts by initializing a two-dimensional matrix of zeros with dimensions 
nbins (the number of bins) and 2. Then, it performs a loop i from 1 to nbins.

For each iteration of the loop, it calculates the mean of the values in the second column (beh(:,2)) 
of the beh matrix, where the rank of the values in the first column (drank) equals i.
 The mean value is stored in the first column of the caf matrix at row i.

It also calculates the mean of the values in the first column (beh(:,1)) of the beh matrix, 
where the rank of the values in the first column (drank) equals i. 
The mean value is stored in the second column of the caf matrix at row i.

This process is repeated for each iteration of the loop, until the loop is finished.
%}

% and plot
figure(5), clf
plot(caf(:,1),caf(:,2),'o-','markerfacecolor','r','markersize',10)
xlabel('Average RT')
ylabel('Average accuracy')

%{
The code above plots the mean reaction time (RT) and mean accuracy for 
each bin (group of values) of the discretized data. 

The CAF (Conditional Accuracy Function) is a graph showing the relationship 
between RT and accuracy.

The plot function is used to create the graph 
and x-axis is labeled "Average RT" while y-axis is labeled "Average accuracy".
%}

%% end
