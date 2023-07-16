%% Matlab for Brain and Cognitive Scientists
% Matlab code for Chapter 30
% Mike X Cohen
% 
% This code accompanies the book 
%      "Matlab for Brain and Cognitive Scientists" (MIT Press, 2017).
% Using the code without following the book may lead to confusion, 
% incorrect data analyses, and misinterpretations of results. 
% Mike X Cohen assumes no responsibility for incorrect use of this code. 

%% define timing parameters

srate = 10000; % sampling rate in Hz

% duration of simulation
sim_dur = 1; % in seconds
timevec = 0:1/srate:sim_dur - 1/srate;

% define input stimulus
input = zeros(1,length(timevec));
input(dsearchn(timevec',.3):dsearchn(timevec',.7)) = 3;
%{
In this section, the timing parameters are defined, including the sampling rate (srate) and 
the duration of the simulation (sim_dur). The timevec variable is created as a time vector 
with the specified duration and sampling rate.

An input stimulus is defined as a vector called input. In this example, a rectangular pulse 
of value 3 is added to the input vector between time points 0.3 and 0.7 seconds.
%}
%% define neuron properties
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
%% run the simulation

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
neuronV(neuronV==volt_reset) = 40;

figure(1), clf
subplot(211)
plot(timevec,neuronV,'k','linew',2)
set(gca,'ylim',[-100 50])
xlabel('Time (s)'), ylabel('Voltage (mV)')

subplot(212)
plot(timevec,input)
set(gca,'ylim',[min(input)-.1 max(input)*1.1])
xlabel('Time (s)'), ylabel('Input current')

%% and now for a population...
%{
Before starting with the code, letfs think about how to do this. Itfs easy
to create 100 neurons; that just involves making a 100-by-time matrix of
membrane potentials instead of the 1-by-time vector that we made for
Reza. But these would all be independent neurons. They need to talk to
each other. In the brain, this is done through chemical and electrical synapses;
in our simulation, we will set the input to each neuron to be the
spikes of other neurons (as in figure 30.1). Spikes from excitatory cells will
bring each neuron closer to its threshold, and spikes from inhibitory cells
will bring each neuron further away from its threshold. 

Because networks are more complicated than individual neurons, itfs
best to take this one step at a time. The first step will be to simulate 100
identical unconnected neurons. We should then see exactly the same
results 100 times (there is no noise in this simulation). This provides a
sanity check: If something goes wrong with 100 independent neurons, we
definitely cannot interpret any results when those neurons are wired
together.
%}
%(re)define timing parameters

srate = 10000; % sampling rate in Hz

N_exc = 80; % number of excitatory neurons
N_inh = 20; % number of inhibitory neurons

% duration of simulation
sim_dur = 1; % in seconds
timevec = 0:1/srate:sim_dur - 1/srate;

% define input stimulus

input = zeros(1,length(timevec));
input(dsearchn(timevec',.3):dsearchn(timevec',.7)) = 3;

% define neuron properties

% voltage parameters
volt_rest   = -70; % resting membrane potential (mV)
volt_thresh = -50; % action potential threshold (mV)
volt_reset  = -75; % post-spike reset voltage (a bit below resting)

% membrane parameters
R_m = 10; % neuron membrane resistance (MOhm)
tau = 10; % time constant of decay (ms)

% initialize neuron membrane potential
neuronV = volt_rest + zeros(N_exc+N_inh,length(timevec));
spiketimes = [];
%{
The neuronV variable is initialized as a vector of the same size as timevec, 
representing the neuron's membrane potential. The spiketimes variable is an 
empty array that will store the time points when the neuron spikes.
%}

% run the simulation

for timei = 1:length(timevec)-1
    
    % test whether spike threshold was passed
    
    spikedNeurons = neuronV(:,timei) > volt_thresh; % threshold exceeded!
    
    % reset voltage
    neuronV(spikedNeurons,timei) = volt_reset;
    %spiketimes = cat(1,spiketimes,timei);
    
    % update voltage in two steps
    %{
    Assuming your code passed the sanity check, itfs time to connect these
neurons into a network. We can use this network as an opportunity to test
how well activity propagates through the network. We will program the
simulation to apply external input only to half of the neurons. To accomplish
this, the input needs to change from a scalar that is applied equally to
all neurons to a vector that allows individualized inputs.
    %}
    % step 1: resting potential plus input scaled by membrane resistance
    restInput = volt_rest + input(timei)*R_m;
    restInput = [restInput*ones(1,(round(N_exc+N_inh)/2)) volt_rest*ones(1,(round(N_exc+N_inh)/2))]';
    
    % step 1.5: integrate output of other neurons
    restInput = restInput + 25*sum(spikedNeurons(1:N_exc)) - 10*sum(spikedNeurons(N_exc+1:end));
    
    % step 2: update voltage with peak
    neuronV(:,timei+1) = restInput + (neuronV(:,timei) - restInput) * exp(-1000/srate/tau);
end
%{
The main simulation loop has been modified to simulate a population of neurons. 
At each time step, the code checks if the membrane potentials of any neurons exceed 
the action potential threshold. If a neuron spikes, its membrane potential is reset to the post-spike reset voltage.

The voltage update is performed in two steps. Step 1 calculates the resting potential plus input scaled 
by the membrane resistance for each neuron. Additionally, the output from other neurons is integrated 
into the input for each neuron, where excitatory neurons have a positive influence (25*sum(spikedNeurons(1:N_exc))) 
and inhibitory neurons have a negative influence (-10*sum(spikedNeurons(N_exc+1:end))).

Step 2 updates the voltage of each neuron using an exponential decay function, similar to the previous example.
%}

% and plot

figure(2), clf

% an image of the membrane potentials of all neurons
subplot(211)
imagesc(timevec,[],-(neuronV-mean(neuronV(:))))
xlabel('Time (s)'), ylabel('Neuron number')
colorbar

% plotting niceties
neuronV(neuronV==volt_reset) = 40;

subplot(212)
h = plot(timevec,neuronV([4 94],:));
xlabel('Time (s)'), ylabel('Voltage (mV)')
set(gca,'ylim',[-80 50])

% question: why do the colors in the top plot change 
% when you re-run this cell?

%% Izhikevich neurons 
% code adapted from Izhikevich 2003
%{
Based on the model proposed by Izhikevich in 2003, this code allows you to simulate 
the behavior of Izhikevich neurons and observe their membrane potential dynamics in 
response to specific input patterns.
%}
% parameters that control neuron's behavior
a = .03; % time scale of recovery variable u. Larger the value of a, quicker the recovery.
b = .25; % sensitivity of the recovery variable u to the sub-threshold fluctuations of the membrane potential v.
c = -60; % after spike reset value of v (caused by fast K+ conductances). Typically -65 mV.
d =   4; % after spike reset value of u (caused by slow Na+ and K+ conductances). Typically value of 2.
%{
The parameters a, b, c, and d control the behavior of the neuron. These parameters 
determine the neuron's intrinsic properties and govern its dynamics. The specific 
values used in this code are chosen to produce a specific type of neuron behavior.
%}
tau = .25; % what is this value in Hz?
tspan = 0:tau:1000;

% define time series of input
T1 = zeros(size(tspan));
T1(dsearchn(tspan',200):dsearchn(tspan',800)) = 1;
%{
The code then defines a time series of input T1. This input can be thought of as an 
external stimulus or input current applied to the neuron. In this example, the input is a 
step function that is active during a specific time interval.
%}

V = -70; % membrane potential
u = b*V; % membrane recovery variable
[VV,uu] = deal(zeros(size(T1)));
%{
The membrane potential (V) and the recovery variable (u) are initialized with their resting values. 
The code then iterates over each time step in tspan and updates the membrane potential and recovery 
variable according to the Izhikevich neuron equations. If the membrane potential reaches a certain 
threshold (V > 30), a spike is generated, and the membrane potential and recovery variable are reset.
%}

for ti=1:length(tspan)
    
    % membrane potential
    V = V + tau*(.04*V^2 + 5*V + 140 - u + T1(ti) ); % T1 is the input
    u = u + tau*a*(b*V-u);
    
    if V > 30 % there was a spike
        VV(ti+1)=30;
        V = c;
        u = u + d;
    else % there was no spike
        VV(ti+1)=V;
    end
    uu(ti+1)=u;
end
%{
In the code snippet you provided, a for loop is used to iterate over each time step (ti) in the 
tspan vector. Within each iteration, the membrane potential (V) and the recovery variable (u) 
are updated according to the Izhikevich neuron equations.

The equation V = V + tau*(.04*V^2 + 5*V + 140 - u + T1(ti)) calculates the change in the 
membrane potential at each time step. It takes into account the current membrane potential (V), 
the recovery variable (u), and the input stimulus (T1(ti)). This equation represents the dynamics 
of the membrane potential based on the specific Izhikevich neuron model.

The equation u = u + tau*a*(b*V-u) updates the recovery variable (u) based on the 
membrane potential (V). This equation accounts for the interaction between the 
membrane potential and the recovery variable, contributing to the overall behavior 
of the neuron.

If the membrane potential exceeds a threshold of 30 mV (V > 30), it is considered a 
spike. In this case, the membrane potential is set to 30 mV (VV(ti+1)=30), representing 
the spike, and the membrane potential and recovery variable are reset to their respective 
values (V = c and u = u + d). This mimics the behavior of an action potential in a neuron.

If the membrane potential does not exceed the threshold (V <= 30), it is considered a 
subthreshold value, and the membrane potential is stored in the VV vector (VV(ti+1)=V). 
This allows the tracking of the membrane potential over time.

The uu vector is used to store the values of the recovery variable (u) at each time step, 
providing additional information about the neuron's dynamics.

Overall, this code performs a simulation of an Izhikevich neuron's behavior over time, 
updating the membrane potential and recovery variable at each time step based on the 
specific equations and conditions of the Izhikevich neuron model.
%}

% and plot it
figure(3), clf
plot(tspan,VV(1:end-1),tspan,10*T1-88);
axis([0 max(tspan) -90 35])
xlabel('Time (ms)'), ylabel('Membrane potential')

%% Rescorla-Wagner-esque learning model
%{
The code snippet provided implements a Rescorla-Wagner-esque learning model. 
This type of model is commonly used to simulate classical conditioning and associative learning.
%}
nTrials = 100;      % number of trials
lrate = .3;         % learning rate

rewProbs = [.7 .2];             % vector representing reward probabilities 

w = .5+zeros(nTrials+1,2);      % weights for each action at each trial (initialized at 0.5).
[action,rewpred,pPickAct1] = deal( zeros(1,nTrials) ); 
% Assign (1, 100) zeros to the following:
% action: vector that stores selected actions for each trial
% rewpred: is a vector that stores the reward prediction errors for each trial.
% pPickAct1: is a vector that stores the probability of picking action 1 for each trial.

for triali=1:nTrials
    % compute probability of picking action
    pPickAct1(triali) = exp(w(triali,1)) / sum(exp(w(triali,:)));
    
    % pick an action based on weighted probability
    action(triali) = 1 + (pPickAct1(triali)<rand);
    
    % is this action rewarded? (convert to -1/+1)
    reward = rand < rewProbs(action(triali));
    
    % compute prediction error (aka delta)
    rewpred(triali) = reward - w(triali,action(triali));
    
    % update weights for the next trial
    w(triali+1,action(triali)) = w(triali,action(triali)) + lrate*rewpred(triali);
    w(triali+1,3-action(triali)) = w(triali,3-action(triali));
end

figure(4), clf
subplot(311), plot(w,'linew',2), legend({'w1';'w2'})
xlabel('Trial'), ylabel('Action weights')
set(gca,'xlim',[1 nTrials])

subplot(312), plot(1:nTrials,rewpred,'linew',2)
hold on, plot(get(gca,'xlim'),[0 0],'k:','linew',3)
xlabel('Trial'), ylabel('Reward prediction error')
set(gca,'xlim',[1 nTrials])

subplot(313), plot(1:nTrials,pPickAct1,'linew',2), set(gca,'ylim',[0 1])
hold on, plot(get(gca,'xlim'),[.5 .5],'k:','linew',3)
xlabel('Trial'), ylabel('Action probability')
set(gca,'xlim',[1 nTrials])

%% end
