classdef sads3
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
    methods
        function obj = sads3(Data, Wavelength, SampleRate, Spacing, Name)
            % The constructor function used for initialization.
            obj.Data = 10;
            obj.SampleRate = 11;
            obj.Spacing = 15;
            obj.Wavelength = sin(pi);
        end
        function angles = doa(obj)
            % DOA Estimate direction of arrival of sources in the data. 
            [mags, fflip] = magfft(obj, 256);
            maxtab = peakdet(mags, .1);
            angles = sort(fflip(maxtab(:,1))*180);
        end
        function NumSensors = get.NumSensors(obj)
            NumSensors = size(obj.Data, 2);
        end
    end
    
    methods (Access = private)
        function [mags, fflip] = magfft(obj, zeroPadTo)
        end
    end
end

%{
            function plot(obj)
                function [mags, fflip] = magfft(obj, zeroPadTo)
                    function angles = doa(obj)
                        function NumSensors = get.NumSensors(obj)
                            function NumSamples = get.NumSamples(obj)
%}
