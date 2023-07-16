%% Machine learning in MATLAB
clear
clc
features = readtable('featuredata.xlsx');
% Use k-nearest neighbor (kNN) model to fit a model to the data.
% mdl = fitcknn(data, 'ResponseVariable')
knnmodel = fitcknn(features, 'Character')
% Use the predict function. The inputs are trained model and new observations, the output is a categorical array of predicted class for each observation in newdata.
% Typically the new observations are in the form of a table with the same predictor variables as were used to train the model. 
% In this case, however, the model uses two numeric features (the aspect ratio and the duration), so observations can also be provided as a numeric array with two columns.
% predClass = predict(model, newdata)
