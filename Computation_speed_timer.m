%% Goal: Create a count down timer that counts user input cycles of operations using tic and toc. Plot the results.
% 1. Long format
format long;
%2. setup the timer and user input
ii = 1;
t = [];
f = input('How many cycles do you want to count? ') +1;
%% 3. Create a condition for the count down
while ii < f 
    tic
    % Generate a vector with a single column that records the time with
    % each loop cycle:
    t(ii, 1) = toc;
    ii = ii + 1;
end
%% 4. Plot the figure
figure
plot(t)
xlabel('number of cycles')
ylabel('computation time')
title('MATLAB Computation Tester')