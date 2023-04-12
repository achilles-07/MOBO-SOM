function output_matrix = som_reduce_dim(data)
% SOM_REDUCE_DIM Reduces the dimension of input matrix from N-D to 2-D using SOM
%   OUTPUT_MATRIX = SOM_REDUCE_DIM(DATA) reduces the dimension of the
%   input matrix DATA from N-D to 2-D using a Self-Organizing Map (SOM) with
%   predefined parameters. The output matrix is a 2-by-M matrix, where M is
%   the number of data points in DATA.

% Set SOM parameters
mapSize = [10 10]; % Size of SOM map (2D)
epochs = 100; % Number of training epochs
initLearnRate = 0.1; % Initial learning rate
finalLearnRate = 0.01; % Final learning rate
initNeighbor = 4; % Initial neighborhood distance
finalNeighbor = 1; % Final neighborhood distance

% Convert input data to a matrix
data = cell2mat(data);

% Create and train SOM
net = selforgmap(mapSize);
net.trainParam.epochs = epochs;
net.trainParam.lr = initLearnRate;
net.trainParam.lr_inc = (finalLearnRate - initLearnRate) / epochs;
net.trainParam.showWindow = false;
net.trainParam.sigma = initNeighbor;
net.trainParam.sigma_dec = (finalNeighbor - initNeighbor) / epochs;
net = train(net,data');

% Compute 2D representation of data
output = net(data');

% Convert output matrix to 2-by-M format
output_matrix = [output(1,:); output(2,:)];
end
