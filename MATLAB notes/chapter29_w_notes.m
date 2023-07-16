%% Matlab for Brain and Cognitive Scientists
% Matlab code for Chapter 29
% Mike X Cohen
% 
% This code accompanies the book 
%      "Matlab for Brain and Cognitive Scientists" (MIT Press, 2017).
% Using the code without following the book may lead to confusion, 
% incorrect data analyses, and misinterpretations of results. 
% Mike X Cohen assumes no responsibility for incorrect use of this code. 

%% fitting a simple function
% Nonlinear model fitting using fminsearch. 

funch = @(t) 1*(t).^2 + 0*t + 0; % try changing the parameters (book sets several implicitly to 0)
%{
The function funch is defined using a function handle @(t) and takes a scalar input t. 
It returns the value of the quadratic function at t. The coefficients of the quadratic 
function can be adjusted by changing the parameters in the function definition.
%}
[xval,funcval] = fminsearch(funch, -2);
%{
The command fminsearch(funch, -2) uses the fminsearch function in MATLAB 
to find the minimum of the function funch with an initial guess of -2 for the 
input argument t. The output xval contains the value of t at which the minimum 
is found, and funcval contains the value of the function funch at the 
minimum value of t.
%}
figure(1), clf
x = -4:.001:4;
plot(x,funch(x),'k','markersize',5), hold on
plot(xval,funcval,'ro','markerface','r'), hold off
axis square
xlabel('X'), ylabel('Y')
%{
Plots the funch function over a range of x values and plots the minimum point 
found by fminsearch on the same graph. The resulting plot should show the 
parabolic shape of the quadratic function and the location of the minimum point.
%}
%% fitting a piece-wise linear model
% This code is fitting a piece-wise linear model to a triangular distribution.
% generate a triangle distribution
a = .2;
c = .6;
b = .9;
x = rand(1,10000);
% The next two lines of code define a piece-wise linear model.
y(x<c) = a + sqrt( x(x<c).*(b-a).*(c-a) );
% The first line calculates the y-values for the left half of the triangle distribution, where x is less than the midpoint c.
y(x>c) = b - sqrt( (1-x(x>c)).*(b-a).*(b-c) );
% The second line calculates the y-values for the right half of the triangle distribution, where x is greater than the midpoint c. 
[y, x] = hist(y,100);
%{
The `hist` function in MATLAB is used to compute and plot a histogram of data. 
In this case, `y` is the data and `x` is the vector of bin centers returned by the function. 
The `hist` function divides the range of `y` into 100 equally spaced bins. The resulting histogram is plotted with 
`x` on the x-axis and the bin counts on the y-axis.
%}
% find a decent starting location
[~, initB] = min(abs(x-0.5));
%{
In the line [~, initB] = min(abs(x-.5));, the function min finds the minimum value of the absolute difference 
between the array x and 0.5, which represents the center of the distribution. The (~) is used to indicate 
that the first output argument of the min function, which is the actual minimum value, is not needed and can be 
discarded. Therefore, the second output argument, which represents the index of the minimum value, is assigned to initB. 
This provides a starting point for the optimization procedure that follows.
%}
% create a function handle
funch = @(initB) fit2segLinear(initB, x, y);
%{
Here, funch is a function handle that takes an input initB and returns the output of the function fit2segLinear(initB,x,y). 
fit2segLinear is a custom function used to fit a piece-wise linear model to the data. It takes the breakpoint initB 
and the x and y data as input arguments. The function fits two separate lines to the data, one for values 
of x less than initB and another for values of x greater than initB. The function then calculates the sum 
of squared errors between the predicted and actual y values for each segment of the piece-wise linear fit. 
This sum of squared errors is the value returned by the function, which is minimized by the fminsearch 
function to find the optimal value of the breakpoint that minimizes the sum of squared errors.
%}
% and fit the model!
figure(2), clf
[optBreakPoint, sse, exitflag, fmininfo] = fminsearch(funch, initB);
%{
This line of code is using the fminsearch function in MATLAB to find the optimal breakpoint in a 
piecewise linear function that best fits a given dataset.

fminsearch is a built-in MATLAB function for finding the minimum of a given function. In this case, 
funch is the function being minimized, which takes the initial breakpoint (initB) as its input argument 
and returns the sum of squared errors (sse) between the actual data and the model predictions based 
on that breakpoint.

The output arguments of fminsearch are:

optBreakPoint: the optimized breakpoint that minimizes the sum of squared errors.
sse: the sum of squared errors for the optimized breakpoint.
exitflag: a flag that indicates the reason for termination of the algorithm. A value of 1 indicates successful convergence.
fmininfo: a structure containing additional information about the optimization process.
%}

%% now with gaussian
%{

In this code, a Gaussian function with added noise is generated and then a nonlinear regression is 
performed to estimate the parameters of the Gaussian function.
%}
peak = 4;
fwhm = 1;
cent = 3;
nois = .5;
%{
The parameters of the Gaussian function are defined at the beginning of the code as 
peak, fwhm, cent, and nois, which correspond to the peak amplitude, full width at 
half maximum, center, and noise level, respectively.
%}
x = -10:0.1:10;
% This creates a row vector x containing evenly spaced values from -10 to 10 with a step size of 0.1.
gaus = peak*exp( -(x-cent).^2 / (2*fwhm^2) );
gaus = gaus + nois*randn(size(gaus));
%{
The Gaussian function is generated using these parameters and stored in the variable gaus. 
Then, some random noise is added to gaus using the randn function and stored back in gaus.
%}

% initialize:  peak  fwhm  center
initParms = [   2     2      -2   ];
funch = @(initParms) fitGaussian(initParms, x, gaus);
%{
Next, the initial parameter values for the Gaussian function are defined as initParms = [2 2 -2], 
which correspond to initial values for peak, fwhm, and center respectively. 
These values are used to initialize the fitting process.

A function handle funch is defined that calls the fitGaussian function, 
passing in the initial parameters as well as the x and gaus data.
%}
figure(3), clf
[outparams, sse, exitflag, fmininfo] = fminsearch(funch, initParms);
%{
Finally, fminsearch is called with funch and initParms as inputs, and the output parameters outparams, 
the sum of squared errors sse, the exit flag exitflag, and the optimization information fmininfo are stored for later use.

outparams, sse, exitflag, and fmininfo are output arguments returned by the built-in function fminsearch. 
They are not user-defined and are specific to the optimization algorithm implemented in fminsearch.
%}
%% show getting stuck in a local minimum using an inverse sinc function
%{
This code defines an inverse sinc function using a function handle and then minimizes the function using fminsearch. 
The minimum value and location of the minimum are stored in xval and funcval, respectively. 
The code then plots the function over a range of x values and highlights the minimum with a red dot.

The purpose of this example is to demonstrate how fminsearch can get stuck in a local minimum instead of 
finding the global minimum. In this case, the inverse sinc function has many local minima, and depending on the 
starting point of the search, fminsearch may converge to a local minimum instead of the global minimum.
%}
% define sinc function and handle
funch = @(x) -sin(x)./x;

% minimize function
[xval, funcval] = fminsearch(funch, 0);

figure(5), clf
x = -50:.01:100;
plot(x, funch(x))
hold on
plot(xval, funcval, 'ro','markersize',10,'markerfacecolor','k')

%% hist function

x = randn(1000, 1);

figure(6), clf
hist(x, 40)

[yy, xx] = hist(x, 40);
% The outputs of hist(x, 40) are assigned to yy and xx, respectively.
hold on
plot(xx, yy, 'r', 'linew', 3)

%% histogram function

x = randn(1000, 1);

figure(7), clf
histogram(x, 40) % this function doesn't exist (yet) in octave

% get outputs and it doesn't plot
hdata = histogram(x, 40);

hold on
xvals = ( hdata.BinEdges(1:end-1) + hdata.BinEdges(2:end) )/2;
%{
hdata.BinEdges gives the edges of the bins as an array. For example, 
if the number of bins is 4, then hdata.BinEdges will be an array with 5 elements, 
indicating the start and end points of each bin.

The expression ( hdata.BinEdges(1:end-1) + hdata.BinEdges(2:end) ) calculates 
the sum of adjacent elements of hdata.BinEdges, resulting in an array of length 
one less than hdata.BinEdges.

Finally, dividing the sum by 2 gives the center point of each bin, 
which is stored in the xvals array.
%}
plot(xvals, hdata.Values, 'r', 'linew', 3)
%{
hdata.Values is a vector that contains the count of observations in each bin of the histogram. 
For example, if hdata.BinWidth is set to 1 and there are 5 observations in the first bin, then hdata.Values(1) will be equal to 5.
%}
%% assigning data points to bins
xidx = zeros(size(x));

for bini=1:hdata.NumBins
    ix = hdata.Data > hdata.BinEdges(bini) & hdata.Data < hdata.BinEdges(bini+1);
    xidx(ix) = bini;
end
%{
The loop iterates over each bin index (bini) in hdata, and creates a logical index ix that is true for 
all data points that fall within the boundaries of that bin. The logical index ix is then used to assign 
the bin index bini to all of the data points that belong to that bin.

After the loop completes, xidx is a vector of the same size as x that contains the bin index for each 
data point in x.
%}
%% end.