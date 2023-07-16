%% Matlab for Brain and Cognitive Scientists
% Matlab code for Chapter 31
% Mike X Cohen
% 
% This code accompanies the book 
%      "Matlab for Brain and Cognitive Scientists" (MIT Press, 2017).
% Using the code without following the book may lead to confusion, 
% incorrect data analyses, and misinterpretations of results. 
% Mike X Cohen assumes no responsibility for incorrect use of this code. 

%% input patterns for backprop learning
%{
Wefll have our model solve the famous gexclusive orh (XOR) problem: Given two inputs
that can independently be on or off, respond TRUE when both are on or
off, and respond FALSE when only one is on (figure 31.1).
%}
% XOR patterns (first two numbers; the third number is for the bias term)
inputs = [ 1 0 1; 0 1 1; 1 1 1; 0 0 1 ];
output = [     0;     0;     1;     1 ];
%{
In the provided code, the inputs variable represents the input patterns,
and the the output vector will act as the teacher to give
the model feedback about its accuracy.
The backpropagation algorithm is applied to train a neural network to correctly classify these 
input patterns and approximate the desired outputs.
%}
% initialization
[outputs, predError] = deal( zeros(1, size(inputs, 1)) );
totalError = nan(1,100);
figure(1), clf

% setup model
nInputsNodes = size(inputs, 2);
nHiddenNodes = 3;
nOutputNodes = 1;
%{
The model setup involves defining the number of nodes in the input layer, hidden layer, and 
output layer of the neural network.
%}

% random initial weights
weights_i2h = randn(nInputsNodes, nHiddenNodes); % random normal distribution with nInputsNodes x nHiddenNodes dimensions. 
weights_h2o = randn(nHiddenNodes, nOutputNodes);
l_rate = 0.3;
%{
The learning rate parameter (l_rate) controls how fast the
weights are adjusted: too low and the model never learns; too high and the
model becomes unstable. It is similar to the learning rate of the Rescorla-
Wagner?like model you saw in the previous chapter. I found, through trialand-
error testing, that a learning rate of 0.2 is pretty good for this model
and for this problem. More sophisticated backpropagation models adjust
their learning rate over time.
%}
% run
toggle = true;
iteration = 0;
max_iterations = 300;

%{
We need two nested loops: one while-loop
in which the model keeps adjusting its weights until it learns or it runs out
of time, and one loop over the four input patterns within each iteration of
the while-loop. The while-loops are generally preferred when you donft know 
the number of iterations (in this case, we stop when learning is finished), 
and for-loops are preferred when you do know the number of iterations 
(in this case, four input patterns). As a precaution to avoid getting stuck 
in an infinite loop, we program in a toggle that breaks out of the while-loop 
if the model never learns the correct input-output mapping. 
Thatfs why we have a variable called max_iter, and exceeding this number 
of iterations will flip the toggle. Then we can run a new simulation.
%}
while toggle    % Continues as long as "toggle" variable is true, and stops when it becomes false.
    
    % loop through inputs
    for ini=1:size(inputs, 1)
        
        % forward part: compute inputs and errors
        
        % input-layer to hidden-layer: raw input -> weights -> sigmoid
        hdLayerResp = inputs(ini,:) * weights_i2h;
        hdLayerResp = 2./(1+exp(-hdLayerResp'*2))-1;
        %{
        The line "hdLayerResp = 2./(1+exp(-hdLayerResp'*2))-1;" applies the activation function to the 
hidden layer response (hdLayerResp). It uses the sigmoid activation function, which transforms any 
number to a value between -1 and +1. The variable ini is the looping index around the four input patterns.
The variable weights_h2o contains the weights from the hidden layer to the output layer.
And weights_i2h contains the weights from the input layer to the hidden layer
        %}
        otLayerResp = hdLayerResp' * weights_h2o;
        otLayerResp = 2./(1+exp(-otLayerResp'*2))-1;
         %{
        Same function as above, but assigns value to "otLayerResp" instead. 
        This completes the forward sweep: the inputs are multiplied by one set
of weights and passed into the hidden layer, and then the hidden layer
activations are multiplied by another set of weights and passed into the
output layer.
        %}
        % prediction error
        predError(ini) = otLayerResp - output(ini);
       
        % collect output-layer responses for plotting
        outputs(ini) = otLayerResp;
        %{
        The line "outputs(ini) = otLayerResp;" stores the output layer response (otLayerResp) for 
the current input pattern in the "outputs" vector. It assigns the value of otLayerResp to the 
corresponding element in the "outputs" vector at index "ini". This allows us to keep track of the 
output values produced by the neural network for each input pattern during the training process.
        %}
        
        % backwards part: adjust weights based on error
        %{
        Now for the second part of the backpropagation algorithm: sending the
error back from the output layer to the hidden layer to the input layer. As
the errors are backpropagated, the weights are adjusted. For each set of
weights, we compute a gdeltah?the amount that the weights should be
changed?and then simply subtract that delta from the weights.
        %}
        % adjust hidden -> output weights
        delta = l_rate * predError(ini) * hdLayerResp;
        %{
        The line "delta = l_rate * predError(ini) * hdLayerResp;" calculates the weight 
adjustment values by multiplying the learning rate (l_rate) with the prediction error 
for the current input pattern (predError(ini)) and the hidden layer response (hdLayerResp).
        %}
        weights_h2o = weights_h2o - delta;
        
        % adjust input -> hidden weights
        backprop = weights_h2o .* (1-hdLayerResp.^2) * inputs(ini,:);
        %{
        Above line calculates the backpropagation term by multiplying the weights between the 
hidden and output layers (weights_h2o) with the derivative of the activation function applied to 
the hidden layer response (1-hdLayerResp.^2), and then further multiplying it by the current input 
pattern (inputs(ini,:)).
        %}
        delta = l_rate * predError(ini) * backprop;
        %{
        The above calculates the weight adjustment values (delta) by multiplying the learning rate (l_rate) 
with the prediction error for the current input pattern (predError(ini)) and the backpropagation term (backprop).
        %}
        weights_i2h = weights_i2h - delta';
        %{
        The above updates the weights_i2h matrix by subtracting the computed delta values. This step applies gradient 
descent to adjust the weights between the input and hidden layers, aiming to minimize the prediction error and improve 
the overall performance of the neural network.
        %}
    end
    
    iteration = iteration+1; % not auto-updated in a while loop!
    totalError(iteration) = sum(predError.^2); % Track total error at each iteration of the loop. 
    
    if totalError(iteration)<.01 || iteration>max_iterations % if total error is less than .01 or if the iteration is greater than set max limit, toggle into false to stop the loop. 
        toggle=false;
    end
    
    % plot
    if mod(iteration, 5) == 0 % checks to see if the current iteration is a multiple of 5. Plot at every 5th iteration of the loop.
        % The mod function has the following format: remainder = mod(dividend, divisor)

        subplot(221), cla
        plot(output, 'bo', 'linew', 2, 'markersize', 20, 'markerfacecolor', 'b')
        hold on
        plot(outputs, 'r*','linew', 2, 'markersize', 10)
        set(gca,'xlim', [0 5], 'ylim', [-.1 1.1], 'xtick', 1:4, 'xticklabel', {'o1';'o2';'o3';'o4'})
        
        subplot(222)
        plot(totalError, 'k', 'linew',2)
        set(gca,'xlim', [0 max_iterations],'ylim', [0 2.4])
        xlabel('Trials'), ylabel('Error^2')
        
        subplot(223)
        imagesc(weights_i2h), set(gca, 'clim', [-2 2])
        title('input -> hidden weights')
        
        subplot(224)
        imagesc(weights_h2o'), set(gca, 'clim', [-2 2])
        title('hidden -> output weights')
        
        drawnow % used to update the figure and display the plots in realtime.
    end % end plotting loop. Note that it is not necessary for plotting to be in a loop. 
    
end % end simulation while-loop

% end
%% k-means clustering
%{
This code performs k-means clustering on a dataset and visualizes the results.
%}
load kmeans_data

% plot the raw data
figure(2), clf
plot(d(:,1), d(:,2),'ko')    % Plot data points from dataset 'd'
axis([0 1 0 1])     % limit x-axis and y-axis in both dimensions to [0, 1]

% k-means clustering
k = 3; % how many clusters?
[groupidx, cents, sumdist, distances] = kmeans(d, k);
%{
The code kmeans(d, k) performs the k-means clustering algorithm on the dataset d with k clusters. 
It assigns each data point in d to one of the k clusters and calculates the centroid (mean) of each cluster.

The output variables are:

groupidx: A vector of size equal to the number of data points, indicating the assigned cluster index for 
each data point. Each element in groupidx represents the cluster assignment for the corresponding data point in d.

cents: A matrix of size k-by-2, where each row represents the coordinates of a centroid for a cluster.

sumdist: The sum of distances squared between the data points and their respective centroids. 

distances: A matrix of size m-by-k, where m is the number of data points and k is the number of clusters. 
Each element (i, j) in distances represents the Euclidean distance between the i-th data point and the j-th centroid.
%}

figure(3), clf

% draw lines from each data point to the centroids of each cluster
lineColors = 'rkb';
hold on
for i=1:length(d)
    plot([ d(i, 1), cents(groupidx(i), 1) ], [ d(i, 2), cents(groupidx(i), 2) ], lineColors(groupidx(i)))
end
%{
The line of code you provided is used to plot a line segment between a data point (d(i, 1), 
d(i, 2)) and the centroid of its assigned cluster (cents(groupidx(i), 1), cents(groupidx(i), 2)). 
The color of the line is determined by lineColors(groupidx(i)), which corresponds to the color 
associated with the cluster assigned to the data point.
%}

% now draw the raw data in different colors
for i=1:k
    plot(d(groupidx==i, 1), d(groupidx==i, 2), [ lineColors(i) 'o' ], 'markerface', 'w')
end
%{
d(groupidx==i, 1) selects the x-coordinates of the data points that belong to cluster i. 
Similarly, d(groupidx==i, 2) selects the y-coordinates of those data points.
[lineColors(i) 'o'] creates a color and marker combination for plotting the data points. 
lineColors(i) retrieves the color associated with cluster i, and 'o' specifies a circle marker.
'markerface', 'w' sets the marker face color to white, creating a hollow marker.
%}
% (why is this code after the previous plotting code?)
plot(cents(:,1), cents(:,2), 'ko', 'markerface', 'g', 'markersize', 10)
set(gca,'xlim',[0 1],'ylim',[0 1])
legend({'group 1';'group 2';'group 3'})
%{
cents(:, 1) retrieves the x-coordinates of the centroids from the cents matrix.
cents(:, 2) retrieves the y-coordinates of the centroids from the cents matrix.
This code is placed after the previous plotting code to add the centroids as a 
final layer on top of the data points and lines. By plotting the centroids separately, 
they are visually distinguished from the data points and clearly indicate the locations 
of the cluster centers in the plot.
%}
%% SVM in-sampling (!!)
%{
 The code implements SVM (Support Vector Machines) for classification using EEG (Electroencephalography) data.
It performs both in-sample testing within the training data and out-of-sample testing from new, unseen data.

In the in-sample testing part, the code trains an SVM model using the EEG data at different time points. 
It then evaluates the trained model by predicting labels using the same training data and computes the 
accuracy. The average accuracy for each time point is calculated and plotted.
%}
% this and the next cell will not work in Octave. try the libsvm toolbox
load EEG_LR
% initialize
accu = zeros(size(timevec));
trueLabels = [ones(size(l_eeg, 3), 1); 2*ones(size(r_eeg, 3), 1)];
%{
ones(size(l_eeg,3), 1) creates a column vector of ones with a size equal to the number of 
trials in the left EEG data. Each one represents a label for the left EEG trials.
2*ones(size(r_eeg,3), 1) creates a column vector of twos with a size equal to 
the number of trials in the right EEG data. Each two represents a label for the right EEG trials.
The semicolon ; concatenates the two column vectors vertically, resulting in a single column 
vector trueLabels that contains the labels for all the EEG trials. The labels for the left EEG trials 
are represented by ones, and the labels for the right EEG trials are represented by twos.
%}
% The for loop iterates over the variable ti starting from 1, incrementing by 20, until it reaches the size of the second dimension of the l_eeg matrix.
for ti = 1:20:size(l_eeg, 2) 
    
    % organize the data for this time point
    data = squeeze( cat(3, l_eeg(:, ti,:), r_eeg(:, ti,:) ))';
    
    % fit the model
    svmModel = fitcsvm(data,trueLabels);
    
    % evaluate the model (note: training data used in testing, 
    % so this is in-sample testing)
    catLabel = predict(svmModel,data);
    
    % average accuracy
    accu(ti) = mean(catLabel==trueLabels);
end

% show a time course of the average accuracy results
figure(4), clf
plot(timevec(1:20:end),accu(1:20:end),'k','linew',2)
set(gca,'ylim',[.3 1])
hold on
plot(get(gca,'xlim'),[.5 .5],'k--')
xlabel('Time (ms)'), ylabel('Accuracy')

% out-sample test
%{
In the out-of-sample testing part, the code performs a more rigorous evaluation. 
It iterates over the time points and trials, removing one trial at a time from the training data and using it as a test trial. 
The SVM model is trained on the remaining trials and then used to predict the label for the test trial. 
The accuracy is computed for each test trial and stored. Finally, the average accuracy over time is calculated and plotted.
%}
trueLabels = [ones(size(l_eeg,3),1); 2*ones(size(r_eeg,3),1)];
accu = zeros(length(timevec),length(trueLabels));

% warning! this takes a while to run...
for ti=1:20:size(l_eeg,2)
    
    % get data for this time point (doesn't change over trials!)
    data = squeeze( cat(3, l_eeg(:,ti,:), r_eeg(:,ti,:) ))';
    %{
    l_eeg(:, ti, :) selects a subset of data from the l_eeg matrix. The : in the third dimension means that 
all elements along that dimension are selected. This subset of data corresponds to the ti-th time point in the l_eeg matrix.

Similarly, r_eeg(:, ti, :) selects a subset of data from the r_eeg matrix 
corresponding to the ti-th time point.

The cat function concatenates the selected subsets of data along the third dimension, 
resulting in a new matrix with dimensions [size(l_eeg, 1), size(l_eeg, 3) + size(r_eeg, 3)]. 
This concatenation combines the information from the l_eeg and r_eeg matrices at the ti-th time point.

The squeeze function removes any singleton dimensions from the concatenated matrix, 
resulting in a new matrix with dimensions [size(l_eeg, 1), size(l_eeg, 3) + size(r_eeg, 3)]. 
This step is performed to ensure that the dimensions of the data matrix are consistent 
and compatible for further processing.

Finally, the transpose operation ' is applied to the concatenated and squeezed matrix, 
resulting in a new matrix with dimensions [size(l_eeg, 3) + size(r_eeg, 3), size(l_eeg, 1)]. 
This transposition swaps the rows and columns of the matrix, making it more convenient 
for subsequent operations.
    %}
    
%{
The loop below performs cross-validation by iteratively training and testing the SVM model on 
different subsets of the data, allowing for the evaluation of the model's performance in 
predicting the labels of unseen test trials.
%}
    for triali=1:length(trueLabels)
    
        templabels = trueLabels;
        
        % remove test trial from training
        traindata = data;
        traindata(triali,:) = [];
        templabels(triali)  = [];
        
        % fit model on training data
        svmModel = fitcsvm(traindata, templabels);
 %{
 The line svmModel = fitcsvm(traindata, templabels) fits a Support Vector
 Machine, where 'fitcsvm' function is used to train an SVM classifier. 
 The SVM model is trained to find an optimal hyperplane that separates the different 
classes in the training data. The traindata matrix contains the features or predictors, 
and the templabels vector contains the corresponding class labels for each training instance
%}
        % Now predict only the trial that was not included
        % in the training set.
        catLabel = predict(svmModel, data(triali,:));
        accu(ti, triali) = catLabel == trueLabels(triali);
    end
end
%{
In the above, predict function is used to make predictions on the test data data(triali,:) 
using the trained SVM model svmModel. The predicted labels are stored in the variable catLabel.

Then, the accuracy of the prediction is computed by comparing the predicted labels (catLabel) 
with the true labels (trueLabels(triali)) using the == operator. If the predicted label matches the 
true label, the accuracy value is set to 1; otherwise, it is set to 0.

The accuracy values are stored in the accu matrix at the corresponding time index ti and trial index triali.
%}

% and plot
figure(4), clf
plot(timevec(1:20:end), mean(accu(1:20:end,:), 2), 'k', 'linew', 6)
%{
accu(1:20:end,:) selects a subset of rows from the accu matrix corresponding to the selected time points. 
This subset contains the accuracy values for each time point and trial.

mean(accu(1:20:end,:), 2) calculates the mean accuracy across trials for each selected time point.
%}
%% end.
