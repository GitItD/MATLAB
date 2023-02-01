% Linear regression of x and y
%% Introduce the data
x = [20.454, 4.68, 6.93, 7.43, 4.58, 6.40, 6.04, 0.846, 3.49, 4.53]';
y = [24.10, 46.8, 69.0, 74.7, 47.5, 63.8, 61.6, 7.76, 0.739, 45.3]';

%% Calculate the slope and graph it
m = x\y;
yA = m*x;
scatter(x, y);
hold on;
plot(x, yA);
grid on;
xlabel('X')
ylabel('Y')

%% Calculate y-intercept
X = [ones(length(x), 1) x];
b = X\y;
yB = X*b;
plot(x, yB, '--');
legend('Scatter X, Y', 'Slope X/Y', 'Slope & Y-Intercept', 'Location', 'best')

%% Find the better coefficient of determination
Slope = 1 - sum((y - yA).^2)/sum((y - mean(y)).^2)
Slope_Yint = 1 - sum((y - yB).^2)/sum((y - mean(y)).^2)

%% Multiple linear regression using least squares
[b, bint] = regress(y, x)