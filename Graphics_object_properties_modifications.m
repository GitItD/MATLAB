% Graphics object properties modifications

%% Make a graph
figure
x = 0:0.01:20;
y = sin(x);
plot(x, y)

%% Set the graphics proerties of the current y and x-axis
figure
hold on
set(gca, 'xcolor', [0.7, 0.7, 0.7])
set(gca, 'ycolor', [0.7, 0.7, 0.7])

%% Redefine "line_f" as the only handle for defining the line of current graphic object
plot(x,y)
line_y = findall(gca, 'type', 'line')
% Address line_f and set it to some new property
set(line_y, 'linewidth', 3, 'linestyle', ':')
% Print the graph exactly as they look (MATLAB attempts to save ink otherwise)
set(gcf, 'InvertHardCopy', 'off')
%% Practice setting a second line on the same graph
z = cos(x);
h = plot(x,z);
prop_h = get(h)
prop_t = get(gcf)
%% Target a specific graphical property and change it
% Set up the figure
figure
hold on
plot(x,y)
plot(x,z);
% Find the blue line and change its properties
find_line_z = findobj(gca, 'color', [0.8500 0.3250 0.0980]);
set(find_line_z, 'color', 'g', 'linewidth', 3);
%% Strange request to turn off the toolbar of the figure
figure
tb = uitoolbar;
tb.Visible = 'off';
%% Create 3 unequal shaped plots on a single figure 
figure
% location = axes('position', [position_x, position_y, length_x, length_y])
h1 = axes('position', [0, 0.8, 1, 0.2])
h2 = axes('position', [0.8, 0, 0.2, 0.8])
h3 = axes('position', [0, 0, 0.8, 0.8])
% Create 3 axes and their handles
set(gcf, 'CurrentAxes', h1)
plot(x,y)
set(gcf, 'CurrentAxes', h2)
plot(x,y)
set(gcf, 'CurrentAxes', h3)
plot(x,y)