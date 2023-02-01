clear
clc

%% Histogram plots
ru = rand(1000, 1);
r = randn(1000, 1);
[y_r, x_r] = hist(r, 50);
[y_ru, x_ru] = hist(ru, 50);
plot(x_r, y_r, 'k', 'linew', 2)
hold on
plot(x_ru, y_ru, 'r', 'linew', 2)
legend({'randn'; 'rand'})
figure
hold on;
histogram(r, 50)
histogram(ru, 50)
%% Subplots
figure
subplot(2, 2, 1)
plot(r)
subplot(2, 2, 2)
plot(ru)
subplot(2, 2, 3)
plot(x_r, y_r)
subplot(2, 2, 4)
plot(x_ru, y_ru)
%% Patch polygons
figure
x = [0 1 1 0];
y = [0 0 1 1];
patch(x, y, 'red')

figure
x2 = [2 5; 2 5; 8 8];
y2 = [4 0; 8 2; 4 0];
patch(x2, y2, 'green')

figure
hold on
patch(x, y, 'red')
patch(x2, y2, 'green')
%% Working with image files
cd(['C:\Users\dun81\OneDrive\Documents\Dun II\ET Online\Library Online\Psychology\Computational Neuroscience\' ...
    'MATLAB for Brain and Cognitive Scientists - The MIT Press (2017)\MATLAB_BrainAndCognitiveScientists-main\ch08'])
figure
pic = imread('saturn.png');
imagesc(pic)
colorchans = {'red'; 'green'; 'blue'};

for ii = 1:3
    subplot(2, 2, ii)       % Create ii number of channels
    imagesc(pic(:,:, ii))           % Create ii number of pages for matrix (see your Excel cheat sheet to understand what this is)
    axis off                    
    set(gca, 'clim', [0 255])               % Create a get current axis and create a color limit ('clim') from 0 to 255
    title([colorchans{ii} ' channel'])          % Create a title which alternates from {'red' = 1}, {'green' = 2}, and {'blue' = 3}
end
subplot(2,2,4), surf(pic(:, :, 1))          % Turn the picture into a 3D surf plot
shading interp                  % Use this to smoothen the figure and allow better visualization
% Running the code below highlights the discinction between contourf and imagesc
figure
subplot(121)
contourf(pic(:,:,1), 20, 'linecolor', 'k')
title('contourf')
subplot(122)
imagesc(pic(:,:,1))
title('imagesc')
% Open the colormap editor to adjust your colors of the image
colormapeditor 

%% Get, set and handle
% Create a tick mark using debug tool to test to see what each code below does.
figure
plot(rand(3));
xlim = get(gca, 'xlim')
yTik = get(gca, 'ytick')
hold on
plot(get(gca, 'xlim'), [0, 0], 'k')

set(gca, 'xlim', [0, 2])
set(gca, 'ytick', [0, 0.5, 0.8, 0.91])
set(gca, 'xlim', [0, 2], 'ytick', [0.5, 1], 'xtick', 0:0.2:2)

line_h = plot(1:10, (1:10).^2);
line_h_return = get(line_h)
gca_return = get(gca)
set(line_h, 'linewidth', 4, 'marker', 'o', ...
    'markeredgecolor', 'k', 'LineStyle', ':', 'Markersize', 25)
plot_hs = zeros(1, 100);
for ii = 1:100
    plot_hs(ii) = plot(randn(max(1,round(rand*10)),1));
end
figure
rand_h = plot(rand(100));
set(rand_h(1:50), 'color', 'k')
set(rand_h(25:75), 'marker', 'o')
set(rand_h([1:10, 20:5:100]), 'linewi', 4)
title('My Art')
gcf_return = get(gcf)
set(gcf, 'color', 'm', 'name', 'My Art')
axes_h = axes;
figure_h = figure;
marker_return = set(rand_h(1), 'marker')

%% Texts in plots
clf
% Put text in designated locaiton on the graph
text(0.6, 0.4, 'Yo!')
for i = 1:1000, text(rand, rand, 'Yo!'), end
text(1.05, 0.7, 'eastside')
% Using graphics handles on texts
clf
txt_h = text(0.5, 0.5, 'Hello World!')
set(txt_h, 'Position', [0.2, 0.7])
set(txt_h, 'color', 'b', 'String', 'Sweden')
% subscripted = underscore = _
% superscripted = caret = ^
% Greek characters = backslash = \
set(txt_h, 'String', 'Z_oidbe^r^g is an \alpha\_crab')
ylabel('Power (\muV^2)')
delete(txt_h)

%% Interacting with plots
figure, plot(randn(100,1)), datacursormode on 
text(0.5, 2.5, 'Click on any point to display data.')
figure, hold on, plot(randn(100,1)), plot(randn(200,1)), plotedit on
text(0.5, 2.5, 'Click on elements of plot to move or delete it.')