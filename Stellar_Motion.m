%% Project 1: Stellar Motion
%{ 
The star spectrum data in the spectra matrix was collected at evenly spaced wavelengths (É…), 
and you are given the starting wavelength (É…start), the spacing (É…delta), and the number of observations.
All values in starData are *1.0e-12
%}
spectra = readtable('starData.xlsx', 'ReadVariableNames', false);
nObs = size(spectra, 1);
lambdaStart = 630.02;
lambdaDelta = 0.14;

%% Task 1
%{
Create a variable named lambdaEnd (É…end) that contains the value of the last wavelength in the recorded spectrum. 
You can calculate lambdaEnd with the expression É…start+(nObs?1)*É…delta. 
Use lambdaEnd to create a vector named lambda (É…) that contains the wavelengths in the spectrum, from É…start to É…end, in steps of É…delta.
%}
lambdaEnd = lambdaStart + (nObs - 1) * lambdaDelta;
lambda = (lambdaStart: lambdaDelta: lambdaEnd);

%% Task 2
%{
Each column of spectra is the spectrum of a different star. The sixth column is the spectrum of star HD 94028. 
Extract the sixth column of spectra to a vector named s.
%}
s = spectra{:,6};

%% Task 3
%{
Plot the spectra (s) as a function of wavelength (lambda). Use point markers (.) and a solid line (-) to connect the points.
Add the x-label "Wavelength" and the y-label "Intensity" to the plot.
%}
plot(lambda, s, 'Marker', '.', 'LineStyle', '-')
xlabel('Wavelength')
ylabel('Intensity')
title('Star Doppler Effect Wavelength and Intensity')

%% Task 4
%{
Recall that the min function allows two outputs, the second of which is the index for the minimum value. 
This index corresponds to the location of the hydrogen-alpha line (É…Ha).
Create two variables, sHa and idx, that contain the minimum value of s and the index of that minimum value.
Find the wavelength of the hydrogen-alpha line by using idx to index into lambda. 
Store the result as lambdaHa (É…Ha).
%}
[sHa, idx] = min(s);
lambdaHa = lambda(idx);

%% Task 5
%{
The point (lambdaHa,sHa) is the location of the hydrogen-alpha line.
Add a point to the existing axes by plotting x = lambdaHa, y = sHa as a red square ("rs") with a marker size ("MarkerSize") of 8.
%}
hold on
plot(lambdaHa, sHa, 'Marker','square','Color','r','MarkerSize',8)

%% Tasks 6
%{
If you zoom in on the plot, you can see that the wavelength of the hydrogen-alpha line of HD 94028 is 656.62 nm, 
which is slightly longer than the laboratory value of 656.28 nm.
Using the hydrogen-alpha wavelength of the star, you can calculate the redshift factor 
(the speed of the star relative to Earth) using the formula z=(É…Ha/656.28)?1. 
You can then calculate the speed by multiplying the redshift factor (z) by the speed of light (299792.458 km/s).

Calculate the redshift factor (z) and the speed (in km/s) at which the star is moving away from Earth. 
Assign the redshift factor to a variable named z and the speed to a variable named speed.
%}
z = (lambdaHa/656.28) - 1;
light_s = 299792.458;
speed = z*(light_s)

%% Project 2: Stellar Motion
[sHa, idx] = min(spectra{:,:});              % We will plot all of spectradata, so we go from [sHa, idx] = min(spectra{:,6}) to just min(spectra{:,:}) which is minimum of all of its data.

% Same calculations as before, but with updated data
lambdaHa = lambda(idx);
z = (lambdaHa/656.28) - 1;
light_s = 299792.458;
speed = z*(light_s)

% Create a for loop that plots the various aspects of the spectra data using different graph properties:
figure
hold on
for ii = 1:7      % Start a for loop ii
    s_plot = spectra{:,ii};                  % Create a variable named "s_plot" that analyzes each column of spectra
    if speed(ii) <= 0                       % And if during the analysis, the speed variable is equal or less than 0...
        plot(lambda, s_plot ,'--')    % Plot "lambda" on x, "s_plot" on y.
    else                                            % Or else...
        plot(lambda, s_plot, 'LineWidth',3)           % Plot them with larger line
    end
end
hold off

legend('HD  30584', 'HD  10032', 'HD  64191', 'HD   5211', 'HD  56030', 'HD  94028', 'SAO102986'); 
xlabel('Wavelength')
ylabel('Intensity')
title('Star Doppler Effect Wavelength and Intensity Comparisons')