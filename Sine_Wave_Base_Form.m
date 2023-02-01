clc
clear

% Play around with the following codes and see how the graph changes.

srate = 1000;   % sampling rate in Hz, this changes how "smooth" the graph curve is
time = 0:1/srate:2;  % units of seconds, this together with frequency determines how many waves per cycle frequency
f = 20;  % units of Hz, this together with time determines how many waves per cycle time
a = 2;  % arbitrary units, this is the amplitude
th = pi/2;  % in radians, phase shift of the position of the wave
sinewave = a*sin(2*pi*f*time+th)   % Formula of sine wave base form
figure
plot(time, sinewave)

% The formula for complex sine waves (csw) is as follows:
csw = exp(i*2*pi*f*time);
figure
% We can plot a 3D graph with a real and an imaginary part:
plot3(time, real(csw), imag(csw))
xlabel('Time (s)')
ylabel('Real part')
zlabel('Imaginary part')
title('A sine wave at 4 Hz')