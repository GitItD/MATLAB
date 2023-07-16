%% Matlab for Brain and Cognitive Scientists
% Matlab code for Chapter 13
% Mike X Cohen
% 
% This code accompanies the book 
%      "Matlab for Brain and Cognitive Scientists" (MIT Press, 2017).
% Using the code without following the book may lead to confusion, 
% incorrect data analyses, and misinterpretations of results. 
% Mike X Cohen assumes no responsibility for incorrect use of this code. 

%% intro to interpolation
%{
A simple example like this is perfect for understanding the behavior of
different algorithms for performing the interpolation. The choice of algorithm
is an optional third input into the function. There are several options,
including linear, three nearest-neighbor variants, and spline.

In practice, you will probably use linear or spline interpolation most
often. The main difference between these is that the spline method is nonlinear
whereas the linear method is linear. When interpolating
over a fairly densely measured grid, the difference between linear and
nonlinear interpolated results will be minor. 
In other words, if there are a lot of data points available, 
the choice of linear or non-linear interpolation may not make 
a big difference in the accuracy of the final result.
%}

% 'measured' data
data = [1 4 3 6 2 19];
datatimes = 1:6;

%{ 
requested (interpolated) time points

We want to up-sample this vector by 1,000 times, as if we were upsampling
from 1 Hz to 1 kHz. The variable newtimes specifies the time
points at which we want to interpolate.  
%}
newtimes  = 1:0.001:6;

%{ 
different interpolation options

This code sets up an array called "interpOptions" that contains 
four different interpolation methods: linear, spline, next and nearest. 
The "for" loop that follows will iterate through this array, so that the
 code inside the loop will be executed once for each interpolation method. 
The "figure(1), clf" line is clearing any previous plots and setting up a new 
figure window for the upcoming plots. The "methodi" variable is used to keep 
track of the current iteration of the loop and will be used to reference the current
 interpolation method from the "interpOptions" array.
%}

interpOptions = {'linear';'spline';'next';'nearest'};

figure(1), clf
for methodi=1:length(interpOptions)
    
    %{ 
define interpolation object
       
    The measurements were taken on a regular lattice (i.e., the spacing between
each successive measurement was the same), so we use griddedInterpolant.

Use griddedInterpolant to perform interpolation on a N-Dimensional gridded data set.
F = griddedInterpolant(datatimes, data, method) creates a 1-D interpolant from a vector of sample points 
"datatimes" and corresponding values "data", with the interpolation method
as the last argument in any syntaxes.     
    %}
    F = griddedInterpolant(datatimes, data, interpOptions{methodi});
    
    %{ 
query that object at requested time points
    
"newdata" is a variable that is assigned the result of the griddedInterpolant function, which takes in three inputs: 
"datatimes", "data", and "interpOptions{methodi}". The function is used to interpolate the data at the specified time 
points in "newtimes" using the specified interpolation method in "interpOptions{methodi}".

The purpose of interpolating the data at these new time points is to increase the sample rate of the data, for example from 1 Hz to 1 kHz. 
The variable "newtimes" is used as an input in the griddedInterpolant function to interpolate the data at the desired time points.
    %}
    newdata = F(newtimes);
    
    %{
    The code below is used for plotting the results of the interpolation 
using different methods specified in the interpOptions variable.
It creates a 2x2 subplot, with one plot for each of the four interpolation methods specified in interpOptions. 
For each subplot, it plots the new interpolated data in black with a thick line, and the original data points as magenta circles. 
It also sets the x-axis limits to [.5 6.5] and sets the title of the plot to the interpolation method being used.

It is logical to sequence the code in this order if the goal is to compare the results of different interpolation methods on the same data.    
    %}
    subplot(2,2,methodi)
    plot(newtimes,newdata,'k','linew',3)
    hold on
    plot(datatimes,data,'go','markersize',15,'markerfacecolor','m')
    set(gca,'xlim',[.5 6.5])
    title([ '''' interpOptions{methodi} '''' ]) % 4 single quotes here to get a single quote of text!
end

%% interpolation vs. extrapolation
%{
This code is an example of using the griddedInterpolant function to interpolate or 
extrapolate data on a 2-D grid. The function takes in the x, y coordinates and the z values 
of the data points, as well as the interpolation method (in this case 'linear' or 'spline'). 
It then creates an object, Flin or Fspl, that can be used to query new data points, in this 
case the missing row 7 of the grid. This way, it compares the different methods of interpolation. 
The code then plots the original data, the interpolated data using the linear method, and the 
interpolated data using the spline method, so you can compare the results.
%}

%{ 
which row to interpolate? Try 7 (interpolation) or 17 (extrapolation)

The variable "point2interp" is used to specify the point at which the data will be interpolated 
or extrapolated. In this case, the value of "point2interp" is set to 17, which means that the 
data will be interpolated or extrapolated at the 17th point.
%}
point2interp = 17;

%{ 
The data in this code is stored in the variable "A" which is the result 
of the built-in Matlab function "peaks" with an input of 10. This function 
creates a 3D surface with peaks, and the input value determines the size of the surface. 
"B" and "C" are also defined as the result of "peaks(10)". 
They are used later in the code for comparison with "A".

The function deal() is used to assign the output of the peaks(10) function to 
multiple variables. Specifically, the output of peaks(10) is assigned to variables A, B, and C. 
The deal() function can be used to assign the elements of an array or the outputs of a function 
to multiple variables in a single line, rather than assigning each output to a variable individually.
%}
[A,B,C] = deal( peaks(10) );

%{ 
define a regular grid but missing row 7. We will then interpolate
or extrapolate using linear and spline methods.

ndgrid is a function in MATLAB that creates a grid of points in N-dimensional space. 
In this case, the function is being used to create a 2-dimensional grid where the x-coordinates 
are defined by the vector [1:6 8:10] and the y-coordinates are defined by the vector 1:10. 
The resulting grid is stored in the variables x and y.

The griddedInterpolant function creates an object that can be used to perform interpolation 
on the data passed as its input. In this case, the input data is the matrix A([1:6, 8:10],:) which is a subset 
of the matrix A created by peaks(10) . The last input 'linear' is the method of interpolation.
The interpolation happens between 1:6, and 8:10, with row 7 missing. The
value of the missing row is being estimated through either linear or spline
interpolatation method.
%}
[x,y] = ndgrid([1:6 8:10],1:10);
Flin  = griddedInterpolant(x,y,A([1:6, 8:10],:),'linear');
Fspl  = griddedInterpolant(x,y,A([1:6, 8:10],:),'spline');

%{ 
query data points

The value of row 7 is estimated using the linear method and stored in the variable B.

The function repmat(point2interp,1,10) creates a matrix where the value of point2interp is 
repeated 1 times in the first dimension (rows) and 10 times in the second dimension (columns). 
This creates a matrix of size [1,10], in which every element is equal to point2interp.
So the resulting matrix is [17 17 17 17 17 17 17 17 17 17].

Meanwhile, (repmat(point2interp,1,10),1:10) is used as the input arguments for the Flin function, 
which is a gridded interpolant of size (1,10) with the value of point2interp. 
In this case, the point2interp is 17, so the function Flin is querying the linear interpolation 
at the point (17,1:10), which is a 1-by-10 vector of interpolated values for the 7th row of the matrix A.

Again, still using linear and spline methods seperately for comparisons. 
%}
B(7,:) = Flin(repmat(point2interp,1,10),1:10);
C(7,:) = Fspl(repmat(point2interp,1,10),1:10);

figure(2), clf
subplot(221)
imagesc(A)
set(gca,'clim',[-7 7])
title('Original')

subplot(222)
imagesc(B)
set(gca,'clim',[-7 7])
title('linear')

subplot(223)
imagesc(C)
set(gca,'clim',[-7 7])
title('spline')

% plot: comparisons of original data (A), linear interpolation (B), and spline interpolatin (C).
subplot(224)
plot(1:10,A(7,:),'k'), hold on
plot(1:10,B(7,:),'b')
plot(1:10,C(7,:),'r')
set(gca,'xlim',[.8 10.2])

legend({'original';'linear';'spline'})

%% an aside on ndgrid
% you'll learn more about making grids and using them for image
% processing in Chapter 26. For now, try to get a feel for what
% the outputs are by running each line and inspecting the outputs.

% no difference...
x = ndgrid(1:3);
x = ndgrid(1:3,2);

% also no difference
x = ndgrid(1:3,[1 2]);
x = ndgrid(1:3,[1 3]);



% now with 2 outputs
[x,y] = ndgrid(1:3);
[x,y] = ndgrid(1:3,2);

% What's the difference between these lines?
[x,y] = ndgrid(1:3,[1 2]);
[x,y] = ndgrid(1:3,[1 3]);


figure(3), clf
[x,y] = ndgrid(1:30,[1:10 21:30]);
subplot(221)
plot(1:20,x), axis square

subplot(223)
imagesc(x), axis square

subplot(222)
plot(y,1:30), axis square

subplot(224)
imagesc(y), axis square

%% 2D interpolation using real EEG data
%{
2D interpolation is used in multichannel neuroscience recordings for topographic 
plotting and to reconstruct activity from electrodes that provided no valid data. 
It helps in visualizing the spatial distribution of activity and facilitates cross-subject data pooling.

Two-dimensional interpolation requires one extra step of complexity,
which is to specify a 2D grid instead of a 1D array of points. To specify a 2D
grid of points, you can use the functions meshgrid or ndgrid.


%}

% load some sample EEG data
load EEGexample.mat

%{ 
convert polar to cartesian coordinates
The variable chanlocs is a structure that contains the channel labels
and locations for each of (in this example) 64 EEG electrodes. It is a format
that was developed for and implemented in the eeglab data analysis toolbox
(Delorme and Makeig 2004). The purpose of the first line of code is to
convert the polar-format electrode locations (radius and angle) to Cartesian-
format locations (x and y coordinates). 

The line of code [eX,eY] = pol2cart(pi/180*[chanlocs.theta],[chanlocs.radius]); is 
converting the polar-format electrode locations (radius and angle) to Cartesian-format 
locations (x and y coordinates) using the pol2cart function. The input to the function is 
the angle in radians (obtained by multiplying the channel location angles in degrees by pi/180) 
and the radius, which are taken from the chanlocs structure that contains the channel labels 
and locations for each of the EEG electrodes. The output of the function is two arrays, eX and eY, 
which contain the x and y coordinates of the electrodes, respectively.
%}
[eX,eY] = pol2cart(pi/180*[chanlocs.theta],[chanlocs.radius]);

%{ 
interpolation factor, and define grid spacing
The next three lines are the points at which we want to specify our grid, which are 100 linearly spaced
points between the minimum and maximum xy positions.
%}
intFact = 100;
interpX = linspace(min(eX),max(eX),intFact);
interpY = linspace(min(eY),max(eY),intFact);

% now define grid in which to interpolate
%{
The function meshgrid supports only 2D and 3D grids and swaps the columns and
rows in the output (this is sometimes useful for image processing, although
it can cause some confusion). The function ndgrid will be used here
because it is applicable in more general situations.

The code creates two 2D matrices, gridX and gridY, using the function ndgrid. 
The function takes two 1D vectors, interpX and interpY, and creates a 2D grid by 
combining all the possible x-coordinates with all the possible y-coordinates.
 The resulting gridX matrix will have the same x-coordinates in each row and the 
gridY matrix will have the same y-coordinates in each column. 
The 2D grid is used to specify the points at which the data is interpolated to 
create a topographic map.

%}
[gridX,gridY] = ndgrid(interpX,interpY);


% shall we have a look at the grids?
figure(4), clf
subplot(222)
imagesc(gridX)

subplot(223)
imagesc(gridY)

subplot(224)
imagesc(gridX+gridY)


%{ 
scatteredInterpolant is a function that creates an interpolant (an object) for a given 
set of scattered data points. The function takes four inputs in this case:
eX' and eY': These are the x and y coordinates of the scattered data points. 
The transpose operator ' is used to ensure that the input is a column vector.

eeg(:,300): This is the data that is being interpolated. In this case, it is a single 
column of data from the eeg matrix at time point 300.

'linear': This is the interpolation method. In this case, linear interpolation is used. 
Other options include 'nearest', 'cubic', etc. 'none': This is the extrapolation method. 
When the requested points fall outside the measurement boundaries, 'none' will set any outside points to NaN (Ågnot a numberÅh).
The scatteredInterpolant function returns an interpolant object F, which can be used to 
evaluate the interpolated data at any point in the grid, using the syntax F(gridX,gridY)

%}
F = scatteredInterpolant(eX',eY',eeg(:,300),'linear','none');
interpDat = F(gridX,gridY);


figure(5), clf

% image the interpolated map
subplot(121)
imagesc(interpX,interpY,interpDat')
axis off, axis image
axis([-.7 .6 -.6 .6])

% draw actual electrode positions
hold on
plot(eX,eY,'ko','markerfacecolor','w','markersize',8)


% this is what the 'map' would look like without interpolation
subplot(122)
scatter(eX,eY,50,eeg(:,300),'filled')
axis image
set(gca,'color','k')
axis([-.7 .6 -.6 .6])

%% using interp1
%{
This code is creating a figure with 4 subplots, each showing the results of using a different interpolation 
method (linear, spline, next, and nearest) to interpolate the given data. 
The original data and the interpolated data are both plotted on each subplot, 
and the x-axis limits are set to .5 and 6.5 to show the difference between the original and interpolated data. 
The code is also adding a legend to the subplots to indicate which data is the original data and which is the interpolated data.

The difference between the method below and the one at the beginning is
that this one uses "interp1" to perform 1-D interpolation on the data and
returns the interpolated data at the desired time points. The second method
uses "griddedInterpolant" to find interpolation on a N-D gridded data set
because the data came from evenly spaced source. 
%}

% Define the 'real' data. Similar to the first cell.
data = [1 4 3 6 2 19];
datatimes = 1:6;
newtimes = 1:0.001:6;

% test four options (linear and spline work in Octave)
interpOptions = {'linear';'spline';'next';'nearest'};

figure(6), clf
for methodi=1:length(interpOptions)
    
    % get the new data
    newdata = interp1(datatimes,data,newtimes,interpOptions{methodi});
    
    % specify subplot
    subplot(2,2,methodi)
    
    % plot data
    plot(newtimes,newdata,'k','linew',3)
    hold on
    plot(datatimes,data,'go','markersize',15,'markerfacecolor','m')
    
    % adjustments..
    set(gca,'xlim',[.5 6.5])
    legend({'interpolated';'original'})
end

%% zero-padding
%{
Because frequency resolution is defined by the number of time points,
the number of time points is also defined by the frequency resolution. I
realize that sounds like circular logic, but it isnÅft, because we can manipulate
the frequency resolution by zero-padding, independent of the time
domain. For example, to double the sampling rate of a time-domain signal
that is N points long, add N zeros to the end of its Fourier spectrum before
taking the inverse Fourier transform. Fortunately, you donÅft need to add
the zeros explicitly; you can specify the N of the inverse Fourier transform
as an optional second input when calling the MATLAB function ifft.
A simple example of applying the zero-padding theorem to perform sincinterpolation
is shown below

The key parameter for the up-sampling is padN. Try changing that
parameter to see the effects on the resulting signal. ItÅfs a good idea to sanity
check the procedure by changing padN to origN, which should have no
effect.
%}

origN =  10;
padN  = 100;

% random data
data = randn(origN,1);

% padded data
datapad = ifft( fft(data)/origN ,padN)*padN;


% and plot!
figure(7), clf
plot(linspace(0,1-1/length(datapad),padN),real(datapad),'bo-','markersize',10)
hold on
plot(linspace(0,1-1/origN,origN),data,'rs-','markersize',10)

legend({'zero-padded';'original'})

%% downsample
%{
Down-sampling data is mostly straightforward (figure 13.8). There is
one additional consideration for down-sampling, which is that higherfrequency
activity can be aliased into lower frequencies.

In many cases, you do not need to apply the low-pass filter yourself.
Many MATLAB functions that perform down-sampling should apply an
anti-aliasing filter prior to removing data. ItÅfs a good idea to check that a
down-sample function you want to use first applies a low-pass filter. If not,
you can apply the low-pass filter yourself before down-sampling.

If you have the MATLAB signal processing toolbox, you can use the
function resample. If not, you can download the (free) Octave package
called Ågsignal.Åh Many other third-party data analysis toolboxes will include
functions to down-sample time series data while applying an anti-aliasing
filter.
%}
origSrate = 1000;
newSrate  =  250;

% specify time vector
time = 0:1/origSrate:1-1/origSrate;

%{ 
smooth the data with a Gaussian

To prevent aliasing from occurring during down-sampling, a low-pass
filter should be applied before removing data points. The cutoff of the filter
should be the Nyquist frequency of the new down-sampled rate (not the
original Nyquist frequency). This is often called an anti-aliasing filter. The
cutoff could also be lower than the new Nyquist frequency, but it should
not be higher.

Most resampling functions do not ask for the sampling rate in hertz or
any other unit. Instead, they take as inputs two integers whose ratio produces
the fraction of the original sampling rate to which you want to downsample.
For example, if the original sampling rate is 1,000 Hz and you want
to down-sample to 250 Hz, you would input the numbers 1 and 4.

Fortunately, you donÅft need to figure out these integers yourself; 
there is a MATLAB function called rat (I'm sure you immediately guessed that r.a.t.
is an acronym for rational fraction approximation). Whenever possible,
down-sample only to a sampling rate that is easily captured by integer fractions.
That is, if the original sampling rate is 1,000 Hz, itÅfs better to downsample
to 250 Hz than to 257.31 Hz. Down-sampling to a sampling rate
that is not an easy integer fraction of the original sampling rate requires
finer interpolation, which increases the risk of inaccuracies.
%}
data = conv( randn(length(time),1), gausswin(400,20) ,'same');

dsFact = origSrate/newSrate;


datads = data(1:dsFact:end);
timeds = time(1:dsFact:end);


figure(8), clf
plot(time,data,'b.-','markersize',10)
hold on
plot(timeds,datads,'ro')

legend({'original';'downsampled'})

%% end
