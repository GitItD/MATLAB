%% MATLAB for Brain and Cognitive Scientists Exercise 29

%% 29.6.2
%{
Look on the Internet for the formula for generating random variables
with a triangular distribution (Wikipedia, for example, should have a
relevant entry). Without looking at the code for this chapter, implement
this formula in MATLAB, in two ways. First, loop through elements
in the vector of random numbers; second, do it without a loop.
Do you also get this existential feeling of wholesome goodness by solving
problems without loops?
%}
% Basic parameters
x = randn(1000, 1);
a = min(x);
b = max(x);
c = mean(x);
% Setup function
f = zeros(size(x));
% Loop through elements of random numbers
for ii = 1:length(x)
    if x(ii)<a || x(ii)>b
        f(ii) = 0;
    elseif x(ii)<=c
        f(ii) = (2*(x(ii)-a))/((b-a)*(c-a));
    else
        f(ii) = (2*(b-x(ii)))/((b-a)*(b-c));
    end
end
% Plot histogram and PDF of f
figure(1), clf
subplot(1,2,1)
histogram(x, 'Normalization', 'pdf')
hold on
plot(x, f, 'r', 'LineWidth', 2)
legend('Histogram', 'Triangular PDF')
% To do it without loop, we set parameters of g. Because there is no if statement,
% we must write the actual inequality properties of the tri distribution.
% Reference <https://mathworld.wolfram.com/TriangularDistribution.html>
g = zeros(size(x));
g(x<a | x>b) = 0;
g(x>a & x<=c) = 2*((x(x>a & x<=c)-a)/((b-a)*(c-a)));
g(x>c & x<b) = 2*((b-x(x>c & x<b))/((b-a)*(b-c)));
% Plot histogram and PDF of g
subplot(1,2,2)
histogram(x, 'Normalization', 'pdf')
hold on
plot(x, g, 'r', 'LineWidth', 2)
legend('Histogram', 'Triangular PDF')
% end
%% 29.6.6
%{
In the triangular-distribution example, perform optimization by Ågbrute
force.Åh This means you should write a loop over all possible x values,
use each x value as a breakpoint and fit the piece-wise linear regression,
and save the R2 for each breakpoint. The best breakpoint parameter is
the one with the highest R2. Confirm that you get the same (or very
similar) answer as with fminsearch. Then run a time test on your
brute-force code and compare that against using fminsearch. This
brute-force method is okay for functions with one parameter. What
would be involved in applying the brute-force method in the Gaussian
example, which has three parameters?
%}
% Create Gaussian parameters
x = randn(1000, 1);
peak = max(x);
fwhm = min(x);
nois = 0.5;
cent = mean(x);
% Perform optimization by writing a loop over all possible values of x
tic
r = zeros(size(x));
g = zeros(size(x));
for ii = 1:length(x);
    if x(ii)==fwhm || x(ii)==peak
        r(ii) = 0;
        g(ii) = 0;
    elseif x(ii)<=cent
        r(ii) = x(ii)<=cent;
    else 
        g(ii) = x(ii)>cent;
    end
end
time_A = toc
% Create Gaussian
tic
gaus = peak*exp( -(x-cent).^2 / (2*fwhm^2) );
gaus = gaus + nois*randn(size(gaus));
% Perform fminsearch using parameters
par = [fwhm, cent, peak];
func = @(par) fitGaussian(par, x, gaus);
time_B = toc
% Sanity check
figure(2), clf
subplot(1,2,1)
plot(x, r, 'ro-', 'LineWidth', 2, 'MarkerSize', 2)
hold on
plot(x, g, 'gs-', 'LineWidth', 2, 'MarkerSize', 2)
legend('Red curve', 'Green curve')
xlabel('X')
ylabel('Amplitude')
title('Comparison of Red and Green curves')
subplot(1,2,2)
[outparams, sse, exitflag, fmininfo] = fminsearch(func, par);
plot(x, gaus, 'b-', 'LineWidth', 2)
xlabel('X')
ylabel('Amplitude')
title('Fitted Gaussian curve')
%end
%% 29.6.9
%{
How robust is the Gaussian parameter search to noise? Re-run the
Gaussian fitting procedure over different noise levels, ranging from 0
to 4 (these are the numbers that scale the output of the randn function
when creating the data). Use at least 10 steps of noise. For each level of
noise, compute the squared differences between the parameters you
specified and the parameters returned by fminsearch. Repeat these
procedures 50 times to be able to take the median over noise levels (you
are likely to get some horrible fits, so median is preferred over mean).
Plot the differences as a function of noise level.
%}
% Generate different noise levels for Gaussian
nois = linspace(0, 4, 10);
gaus_noisy = zeros(size(gaus,1), length(nois));
for ii = 1:length(nois)
    gaus_noisy(:,ii) = gaus + nois(:,ii)*randn(size(gaus));
end
%% 29.6.11
%{
Weibull functions are occasionally used in psychophysics for fitting
behavioral data. Find the formula for the Weibull function online and
implement it in MATLAB. Try different parameters until you find a
curve that you find visually pleasing. Next, come up with six data
points that are near that curve (donÅft be too exact?you should be able
to fit the function to the data). Write a function to fit a Weibull function
and return the two parameters based on the six data points that
you simulated.
%}
%% 29.6.15
%{
Here are 6 data points: [0 0.05 0.3 0.7 0.9 1] (these are the y-axis values;
the x-axis values are 1:6); fit these to a sigmoid and determine the
parameters. Then do it again adding +1 to all the data values. Plot the
data and the fitted sigmoid on the same graph. Make sure to select
the starting parameters carefully.
%}