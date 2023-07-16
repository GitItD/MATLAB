classdef sads2
    % Sensor array data set class
    properties (GetAccess = private)
        Wavelength      % Wavelength of source (m)
    end
    properties (Constant)
        c = 3e8;               % Speed of wave in medium (m/s)
    end
    properties (Dependent)
        NumSensors     % Number of sensors
        Numsamples      % Number of samples
    end
    properties
        Data                        %samples sensor data
        Spacing                 % spacing of array (m)
        SampleRate          % Sample rate (Hz)
        Name                    % Sensor array test run name
    end
end