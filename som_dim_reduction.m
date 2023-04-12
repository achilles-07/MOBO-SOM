function [mapped_data, net] = som_dim_reduction(data)

% Specify the dimensions of the SOM grid
x_dim = 10;
y_dim = 10;

% Create the SOM network
net = selforgmap([x_dim, y_dim]);

% Train the SOM network using the input data
net = train(net, data');

% Map the input data to the 2D SOM grid
mapped_data = net(data');

% Plot the 2D SOM grid
figure;
plotsompos(net, mapped_data);

end