%% MATLAB for Brain and Cognitive Scientists Exercise 30

%% 30.5.1
%{
Reza was given DC input. Adjust the code that creates figure 30.2 to
provide oscillatory input. Test several different frequencies. Given the
amount of time that is simulated, what is the slowest frequency that is
sensible to use?
%}

% define timing parameters

srate = 10000; % sampling rate in Hz

% duration of simulation
sim_dur = 1; % in seconds
timevec = 0:1/srate:sim_dur - 1/srate;

% define input stimulus
input = zeros(1,length(timevec));
input(dsearchn(timevec', 0.3):dsearchn(timevec',.7)) = 3;
%{
In this section, the timing parameters are defined, including the sampling rate (srate) and 
the duration of the simulation (sim_dur). The timevec variable is created as a time vector 
with the specified duration and sampling rate.

An input stimulus is defined as a vector called input. In this example, a rectangular pulse 
of value 3 is added to the input vector between time points 0.3 and 0.7 seconds.
%}

% define neuron properties
%{
We will create a model neuron and give it some input. I think Reza is a
good name for our simulated IAF neuron. We start by initializing Rezafs
parameters. Reza will have his resting, spike threshold, and reset values set
at -70, -50, and -75 millivolts (mV), respectively. The membrane resistance
is 10 megaohms (Mƒ¶), and the decay time constant is 10 milliseconds
%}

% voltage parameters
volt_rest   = -70; % resting membrane potential (mV)
volt_thresh = -50; % action potential threshold (mV)
volt_reset  = -75; % post-spike reset voltage (a bit below resting)

% membrane parameters
R_m = 10; % neuron membrane resistance (MOhm)
tau = 10; % time constant of decay (ms)

% initialize neuron membrane potential
neuronV = volt_rest + zeros(size(timevec));
spiketimes = [];
%{
In this section, the properties of the IAF neuron are defined. 
This includes the voltage parameters such as the resting potential (volt_rest), 
action potential threshold (volt_thresh), and post-spike reset voltage (volt_reset). 
The membrane parameters are also specified, including the membrane resistance 
(R_m) and the time constant of decay (tau).

The neuronV variable is initialized as a vector of the same size as timevec, representing 
the neuron's membrane potential. The spiketimes variable is an empty array that will store 
the time points when the neuron spikes.
%}

% run the simulation

for timei = 1:length(timevec)-1
    % test whether spike threshold was passed
    
    if neuronV(timei) > volt_thresh % threshold exceeded!
        % reset voltage
        neuronV(timei) = volt_reset;
        spiketimes = cat(1,spiketimes,timei);
    end
    
   % update voltage in two steps
    
    % step 1: resting potential plus input scaled by membrane resistance
    restInput = volt_rest + input(timei)*R_m;
    
    % step 2: update voltage with peak
    neuronV(timei+1) = restInput + (neuronV(timei) - restInput) * exp(-1000/srate/tau);
    
end
%{
In the main simulation loop, the membrane potential of the neuron is updated 
at each time step. First, the code checks if the neuron's membrane potential exceeds 
the action potential threshold. If it does, the membrane potential is reset to the post-spike 
reset voltage, and the spike time is recorded in the spiketimes array.

Next, the neuron's membrane potential is updated in two steps. In step 1, the resting 
potential is added to the input scaled by the membrane resistance. In step 2, the 
membrane potential is updated using exponential decay based on the time constant (tau).
%}

% plotting niceties
neuronV(neuronV == volt_reset) = 40;

figure(1), clf
plot(timevec, neuronV, 'k', 'linew', 2)
set(gca, 'ylim', [-100 50])
xlabel('Time (s)'), ylabel('Voltage (mV)')
% End
%% 30.5.8
%{
Using the Izhikevich code for a single neuron, test how the neuron
responds to oscillatory inputs of different amplitudes. Using two loops
(three including the loop over simulation time), stimulate the neuron
using sine waves that vary from 1 Hz to 60 Hz and amplitudes that vary
from 1 to 30. Youfll need to add an offset to the sine wave input to
avoid negative values. For each 1-second simulation, count the number
of action potentials and show the results in a 2D plot of stimulation
frequency by amplitude. Do the results look different when using different
neuron parameters?

Tasks:
1. Change amplitudes
2. 3 loops, simulate neuron w/ sine that vary from 1 - 60 Hz, amplitude 1 - 30
3. Count action potential/1 sec
%}
% parameters that control neuron's behavior
a = 0.03; % time scale of recovery variable u. Larger the value of a, quicker the recovery.
b = 0.2; % sensitivity of the recovery variable u to the sub-threshold fluctuations of the membrane potential v.
c = -65; % after spike reset value of v (caused by fast K+ conductances). Typically -65 mV.
d =   2; % after spike reset value of u (caused by slow Na+ and K+ conductances). Typically value of 2.

% parameters for sin wave
freq = 1:60;   % Frequencies to test
amp = 1:30;    % Amplitudes to test

% parameters for time and setup spikecounter
tau = 0.25; % what is this value in Hz?
tspan = 0:tau:1000;
spikecount = zeros(length(amp), length(freq));

for ai = 1:length(amp)
    amplitude = amp(ai);
    
    for fi = 1:length(freq)
        frequency = freq(fi);
        
        % Generate the input sine wave with offset
        input = amplitude * sin(2*pi*frequency*(tspan/1000)) + 10;
        % Initial voltage
        V = -70;    % membrane potential
        u = b*V;    % membrane recovery variable
        for ti = 2:length(tspan)    % Offset to avoid negative values
            % membrane potential
            V = V + (0.04 * V^2 + 5 * V + 140 - u + input(ti-1)) * tau;
            u = u + a*(b*V-u);
            if V >= 30 % there was a spike
                VV(ti+1) = 30;
                V = c;
                u = u + d;
                spikecount(ai, fi) = spikecount(ai, fi) + 1;
            else % there was no spike
                VV(ti+1) = V;
            end
        end
    end
end
    
    % and plot it
    figure (1),clf
    contourf(freq, amp, spikecount, 40, 'linecolor','none')
    colormap hot, colorbar;
    xlabel('Frequency (Hz)'), ylabel('Amplitude')
    % end
    %% 30.5.8 ChatGPT
    % Neuron parameters
    a = 0.02;
    b = 0.2;
    c = -65;
    d = 8;
    
    % Simulation parameters
    dt = 0.1;      % Time step
    T = 1000;      % Total simulation time in milliseconds
    t = 0:dt:T;    % Time vector
    N = numel(t);  % Number of time steps
    
    % Input parameters
    frequencies = 1:60;   % Frequencies to test
    amplitudes = 1:30;    % Amplitudes to test
    numFreq = length(frequencies);
    numAmp = length(amplitudes);
    spikeCount = zeros(numAmp, numFreq);  % Spike count for each amplitude and frequency combination
    
    % Loop over amplitudes and frequencies
    for iAmp = 1:numAmp
        for iFreq = 1:numFreq
            frequency = frequencies(iFreq);
            amplitude = amplitudes(iAmp);
            
            % Generate the input sine wave with offset
            input = amplitude * sin(2*pi*frequency*t/1000) + 10;
            
            % Initialize neuron variables
            v = -70 * ones(N, 1);  % Initial membrane potential
            u = b * v;             % Initial recovery variable
            spikeTimes = [];       % Spike times
            
            % Simulate the neuron dynamics
            for i = 2:N
                dv = (0.04 * v(i-1)^2 + 5 * v(i-1) + 140 - u(i-1) + input(i-1)) * dt;
                du = (a * (b * v(i-1) - u(i-1))) * dt;
                
                v(i) = v(i-1) + dv;
                u(i) = u(i-1) + du;
                
                if v(i) >= 30
                    v(i) = c;
                    u(i) = u(i) + d;
                    spikeTimes = [spikeTimes, t(i)];  % Store spike times
                end
            end
            
            % Count the number of spikes in 1 second
            spikeCount(iAmp, iFreq) = sum(spikeTimes >= 900 & spikeTimes < 1000);
        end
    end
    
    % Plot the results
    figure(2), clf
    imagesc(frequencies, amplitudes, spikeCount);
    xlabel('Frequency (Hz)');
    ylabel('Amplitude');
    title('Spike Count');
    colorbar;