classdef sads
    % Sensor array data set class
    properties 
        Wavelength      % Wavelength of source (m)
        c = 3e8;               % Speed of wave in medium (m/s)
        NumSensors     % Number of sensors
        Numsamples      % Number os samples
        Data                        %samples sensor data
        Spacing                 % spacing of array (m)
        SampleRate          % Sample rate (Hz)
        Name                    % Sensor array test run name
    end
end