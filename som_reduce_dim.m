function output = som_reduce_dim(data)

% Set SOM parameters
mapSize = [10 10]; % Size of SOM map (2D)
epochs = 100; % Number of training epochs
initLearnRate = 0.1; % Initial learning rate
finalLearnRate = 0.01; % Final learning rate
initNeighbor = 4; % Initial neighborhood distance
finalNeighbor = 1; % Final neighborhood distance

% Initialize SOM weights
weights = rand([prod(mapSize), size(data, 2)]);

% Train SOM
for epoch = 1:epochs
    % Update learning rate and neighborhood distance
    learnRate = initLearnRate + (epoch-1) * (finalLearnRate - initLearnRate) / (epochs-1);
    neighbor = initNeighbor + (epoch-1) * (finalNeighbor - initNeighbor) / (epochs-1);

    % Update weights for each data point
    for i = 1:size(data, 1)
        % Find the best-matching unit (BMU)
        distances = sum((weights - data(i,:)).^2, 2);
        [~, bmu] = min(distances);

        % Update weights for BMU and its neighbors
        [row, col] = ind2sub(mapSize, bmu);
        for j = 1:prod(mapSize)
            [r, c] = ind2sub(mapSize, j);
            d = sqrt((r-row)^2 + (c-col)^2);
            if d <= neighbor
                % Update weight using learning rate and neighborhood function
                weights(j,:) = weights(j,:) + learnRate * (data(i,:) - weights(j,:));
            end
        end
    end
end

% Compute 2D representation of data
output = zeros(size(data,1), 2);
for i = 1:size(data, 1)
    % Find the BMU
    distances = sum((weights - data(i,:)).^2, 2);
    [~, bmu] = min(distances);

    % Convert BMU index to 2D coordinates
    [row, col] = ind2sub(mapSize, bmu);

    % Save 2D coordinates to output matrix
    output(i,:) = [col, row];
end

% Plot 2D representation of data
plot(output(:,1), output(:,2), '.');
