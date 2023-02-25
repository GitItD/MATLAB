clear
clc
%% 14.6.1
%{
Combine figures 14.3 and 14.4 to plot dŒ and response bias on the same
hits and false alarms space. Do the lines look like they are orthogonal?
Looking at figure 14.2, it could be that the hit and false alarm rates are
better conceptualized in log space rather than linear space. Adjust the
plots so the axes are logarithmically scaled (hint: set). How does it
look? Perhaps try having only one of the axes in logarithmic scale.
%}

% Get data (hit and miss rate)
x  = .01:.01:.99;
% Get normal curves of data
dp = bsxfun(@minus, norminv(x)', norminv(x));
% Get d' 
dp2 = [1, 1.5, 2, 2.5];
tol = .01;
% Plot d'
figure(1), clf
contourf(x, x, dp, 80, 'linecolor', 'none')
xlabel('FPR')
ylabel('TPR')
title('Discrimination and Response Bias in TPR and FPR')
hold on
axis square
col = 'rbmk';

for di = 1:length(dp2)
    id = find(dp>dp2(di)-tol & dp<dp2(di) + tol);
    [yi, xi] = ind2sub(size(dp), id);
    plot(x(xi), x(yi), [col(di) 'o-'], 'linew', 4, 'markersize', 9)
end

% Get normal rb
rb = -bsxfun(@plus, norminv(x)', norminv(x))/2;
% Get rb
rb2 = [.3, .5, .9, 1.5];
% Plot dp and rb
for ri = 1:length(rb2)   
    ir = find(rb>rb2(ri) - tol & rb<rb2(ri) + tol);
    [ai, bi] = ind2sub(size(rb), ir);
    plot(x(ai), x(bi), [col(ri) 'o-'], 'linew', 4, 'markersize', 9)
end

% Adjust the plots so the axes are log scaled (hint: set).
figure(2), clf
contourf(x, x, dp, 80, 'linecolor', 'none')
xlabel('FPR')
ylabel('TPR')
title('Discrimination and Response Rate in TPR and FPR Semi-log Scaled')
hold on
set(gca, 'YScale', 'log')
col = 'rbmk';

for di = 1:length(dp2)
    id = find(dp>dp2(di)-tol & dp<dp2(di) + tol);
    [yi, xi] = ind2sub(size(dp), id);
    plot(x(xi), x(yi), [col(di) 'o-'], 'linew', 4, 'markersize', 9)
end
for ri = 1:length(rb2)
    ir = find(rb>rb2(ri) - tol & rb<rb2(ri) + tol);
    [ai, bi] = ind2sub(size(rb), ir);
    plot(x(ai), x(bi), [col(ri) 'o-'], 'linew', 4, 'markersize', 9)
end

%% 14.6.2
%{
Repeat the 2D dŒ plots, but use imagesc instead of contourf. How do
the results change, and do you need to make any adjustments to make
the results accurate?
%}

% Get data (hit and miss rate)
x  = .01:.01:.99;
% Get normal curves of data
dp = -bsxfun(@plus, norminv(x)', norminv(x));
% Get d' 
dp2 = [1, 1.5, 2, 2.5];
tol = .01;
% Plot d'
figure(3), clf
imagesc(x, x, dp)
xlabel('FPR')
ylabel('TPR')
hold on
axis square

for di = 1:length(dp2)
    id = find(dp>dp2(di)-tol & dp<dp2(di) + tol);
    [yi, xi] = ind2sub(size(dp), id);
    plot(x(xi), x(yi), [col(di) 'o-'], 'linew', 4, 'markersize', 9)
end

%% 14.6.3
%{
Now reproduce exercises 1 and 2, but initialize the x-values to be z-units
instead of proportion units. If you make the adjustment correctly, you
should get straight diagonal lines instead of parabolas. (For bonus
points, do this by deleting two norminvfs and adding one).
%}

% Get data (hit and miss rate)
x  = .001:.003:.99;
x = norminv(x, 0, 1);       % The definition of x-values to z-units is to have mean = 0 and standard deviation = 1
% Get normal curves of data
dp = bsxfun(@minus, (x)', (x));
% Get d' 
dp2 = [1, 1.5, 2, 2.5];
tol = .01;
% Plot d'
figure(4), clf
imagesc(x, x, dp)
xlabel('FPR')
ylabel('TPR')
title('Initiate x-values into z-units')
hold on
axis square

for di = 1:length(dp2)
    id = find(dp>dp2(di)-tol & dp<dp2(di) + tol);
    [yi, xi] = ind2sub(size(dp), id);
    plot(x(xi), x(yi), [col(di) 'o-'], 'linew', 4, 'markersize', 4)
end

%% 14.6.5
%{
I used bsxfun in section 14.3 because I donft like loops. Reproduce the
same result using two for-loops. Then use the tic/toc function pair to
check the timing. How much extra time do you need to use for-loops
instead of bsxfun? Now obtain the same result using repmat (no
loops, but no bsxfun). How does the timing of repmat compare to
that of for-loops and of bsxfun?
%}

% Get data (hit and miss rate)
x  = .01:.01:.99;
% Get normal curves of data
tic
dp = -bsxfun(@plus, norminv(x)', norminv(x));
time_bsxfun = toc

tic
for dpi = 1:length(x)
    -(norminv(x(dpi))' + norminv(x(dpi)));
end
time_forloops = toc

%{
tic
dpx = -(norminv(x(dpi))' + norminv(x(dpi)))
repmat(dpx, [1:length(x), 1])
time_repmat = toc
%}

% Get d' 
dp2 = [1, 1.5, 2, 2.5];
tol = .01;
% Plot d'
figure(5), clf
imagesc(x, x, dp)
xlabel('FPR')
ylabel('TPR')
hold on
axis square

for di = 1:length(dp2)
    id = find(dp>dp2(di)-tol & dp<dp2(di) + tol);
    [yi, xi] = ind2sub(size(dp), id);
    plot(x(xi), x(yi), [col(di) 'o-'], 'linew', 4, 'markersize', 9)
end