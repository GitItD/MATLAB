%% MATLAB for Brain and Cognitive Scientists Exercise 30

%% 31.4.2
%{
Run k-means clustering on the 2D data for 2 to 8 groups, 100 times per k.
Devise a method to determine whether the groupings are stable (e.g.,
you might compute the average distances to the nearest vs. all other
centroids). Plot this stability measure as a function of k. This exercise
shows that k-means clustering can be a powerful technique but can
also go awry when suboptimal parameters are used.
%}
% Run the k-means algorithm for k values ranging from 2 to 8, and repeat the process 100 times for each value of k.
% Generate random numbers from a normal distribution
mu = 5; % mean
sigma = 1; % standard deviation
n = 100; % number of random numbers to generate
rnum = sigma * randn(n, 2) + mu; % random number generated in 2 columns
% Scale and shift the random numbers to the desired range of between 2 - 8
lowerBound = 2;
upperBound = 8;
scaledNumbers = (rnum - min(rnum(:))) / (max(rnum(:)) - min(rnum(:)));
d = lowerBound + (upperBound - lowerBound) * scaledNumbers;
% create a 2D figure
figure(1), clf
plot(d(:,1), d(:,2), 'ko')
hold on;
% create k-clusters
k = 3;
[gidx, cent, sumd, dis] = kmeans(d, k);
gscatter(d(:,1), d(:,2), gidx);
plot(cent(:,1), cent(:,2), 'k*', 'markersize', 10);
legend('Cluster 1', 'Cluster 2', 'Cluster 3', 'Centroids');
% end
%% 31.4.3
%{
Can you replace these two lines with one line? You should do this by
extracting N - 1 trials from data, rather than deleting one trial in
traindata. Make sure your code works for any trial, not only the first
or last trial (that would be too easy).

traindata = data;
traindata(triali,:) = [];
%}

% The code takes too long to run, therefore we'll skip this problem.

% end
%% 31.4.7
%{
When XOR is represented as zeros and ones, it can be solved without a
fancy nonlinear backpropagation model, using only one line of code.
Figure out how to reproduce the correct output based on the first two
columns of the input. Of course, this happens to work for this particular
instantiation of XOR; you shouldnÅft think that all problems solvable
with backpropagation can be solved with a single line of code.
%}
% I think we may be able to develop a backpropagation model without
% backpropagation if the predicted error is always 0. To do this, we need
% to make the output always equal to 1. 

p = [0 1 2 3 4 5];
t = [0 0 0 1 1 1];

net = feedforwardnet(2, 'trainrp');

net.trainParam.epochs = 50;
net.trainParam.show = 10;
net.trainParam.goal = 0.1;
net = train(net, p, t);
a = net(p)
%<Source: https://www.mathworks.com/help/deeplearning/ref/trainrp.html>
% end
%% 31.4.7 ChatGPT
% XOR input data
input = [0 0; 0 1; 1 0; 1 1];

% Compute XOR output using logical XOR operator
output = xor(input(:, 1), input(:, 2));

% Plot XOR input and output
figure(2), clf
subplot(1, 2, 1);
scatter(input(:, 1), input(:, 2), 100, 'filled');
xlabel('Input 1');
ylabel('Input 2');
title('XOR Input');
axis([-0.5 1.5 -0.5 1.5]);

subplot(1, 2, 2);
scatter(input(:, 1), input(:, 2), 100, output, 'filled');
xlabel('Input 1');
ylabel('Input 2');
title('XOR Output');
axis([-0.5 1.5 -0.5 1.5]);
colormap([0 0 1; 1 0 0]);
colorbar('Ticks', [0, 1], 'TickLabels', {'0', '1'});
% end
%% 31.4.8
%{
Another option for an activation function is the hyperbolic tangent
(MATLAB function tanh). Your mission here is to use a sigmoid
function (see the formulation in exercise 14 of chapter 29) to obtain
the exact results as the following code. YouÅfll need to figure out how to
set the sigmoid parameters; depending on your level of math background,
you can try this either analytically or empirically.
%}
% The code result you need to obtain:
x = -3:0.1:3;
af = tanh(x);
% The Sigmoid function:
m = (min(x) + max(x)) / 2;   % x-value of function midpoint
s = 1.5;      % steepness of function
v = max(x);     % maximum value of Sigmoid
sig_f = v ./ (1 + exp(-s*(x - m)));
figure(3), clf
hold on
plot(x, af), plot(x, sig_f, 'ro')
legend('af', 'sig_f')
% End
%% 31.4.8 ChatGPT + Dun_Ning input
x = -3:0.1:3;
af = tanh(x);
sig_f = (exp(x) - exp(-x)) ./ (exp(x) + exp(-x));
figure(4), clf
hold on
plot(x, af), plot(x, sig_f, 'ro')
legend('af', 'sig_f')
% End
%% 31.4.9
%{
In the backpropagation model, I guess you played around a bit with the
learning rate. LetÅfs make this more quantitative. Put the model into a
function that takes a learning rate as input and outputs whether the
model successfully learned or stopped because it exceeded its maximum
iteration limit (you might want to turn off plotting). Then run
the model 50 times using 30 different learning rate parameters between
0 and 1 (thus, 1,500 simulations). Plot average successes as a function
of learning rate.
%}
% End