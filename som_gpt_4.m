function [Y, W] = som(data)
% SOM function for dimensionality reduction
% Input:
% data - a matrix of N-dimensional points (each point is a row)
% Output:
% Y - a matrix of 2-D points (each point is a row)
% W - the weight matrix of the trained SOM

% Set SOM parameters
map_size = [10, 10]; % size of the 2-D SOM map
epochs = 500; % number of iterations for SOM training
n = size(data, 2); % dimensionality of input data
m = prod(map_size); % number of neurons in SOM map
alpha_0 = 0.01; % initial learning rate
sigma_0 = max(map_size)/2; % initial neighborhood radius
tau1 = epochs/log(sigma_0); % time constant for neighborhood radius
tau2 = epochs; % time constant for learning rate

% Initialize SOM weight matrix randomly
W = randn(n, m);

% Normalize input data to zero mean and unit variance
data = (data - mean(data)) ./ std(data);

% Train SOM
for t = 1:epochs
    % Pick a random input data point
    i = randi(size(data, 1));
    
    % Compute distances between input data point and SOM neurons
    diff = repmat(data(i,:), m, 1)' - W;
    dist = sum(diff.^2, 1);
    
    % Find the winning neuron (i.e. the one with the smallest distance)
    [~, j] = min(dist);
    
    % Update the weight vectors of the SOM neurons
    [jx, jy] = ind2sub(map_size, j); % convert linear index to 2-D index
    k = 1:m;
    [kx, ky] = ind2sub(map_size, k); % convert linear index to 2-D index
    d = sqrt((jx-kx).^2 + (jy-ky).^2); % distance between neurons
    h = exp(-d.^2/(2*sigma_0^2)); % neighborhood function
    W = W + alpha_0*h.*(diff./(sigma_0^2)); % weight update
    
    % Update learning rate and neighborhood radius
    alpha = alpha_0*exp(-t/tau2);
    sigma = sigma_0*exp(-t/tau1);
end

% Map input data to 2-D points
Y = zeros(size(data, 1), 2);
for i = 1:size(data, 1)
    diff = repmat(data(i,:), m, 1)' - W;
    dist = sum(diff.^2, 1);
    [~, j] = min(dist);
    [jx, jy] = ind2sub(map_size, j);
    Y(i,:) = [jx, jy];
end

end
