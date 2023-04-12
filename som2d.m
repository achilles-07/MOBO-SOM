function y = som2d(Y)
    for i=1:length(Y)
            X(i,:) = Y(i).CostOriginal';
    end
   
%     X = Y;
    [nrows ncols] = size(X);

    %% Initializing Parameters of SOM
    % Weight matrix
    W = rand(nrows,ncols);
    % Learning Rate
    eta0 = 0.5;
    % Standard Deviation
    T = X';
    x = T(1,:);
     y = T(2,:);
    dist = sqrt((x).^2 + (y).^2);
    sigma0 = std(dist);
    % Max number of Iterations
    tau = 50;
    %% SOM Algorithm
    for iter = 1:tau
        eta = eta0 * exp(-iter/tau);
        sigma = sigma0*exp(-iter/tau);
        % Iterating through input
        for i=1:size(X(:,1))
        [n I] = pdist2(W,X(i,:),'euclidean','Smallest',1);
            for j = 1:size(W(:,1))
                h(j,:) = exp((-((norm(W(j,:) - W(I,:))))^2)/(2*sigma*sigma));
%                 if(isnan(h(j,:)))
%                     h(j,:) = 0;
%                 end
            end
            % Update Weights
            W = W + eta.*h.*(X-W);
        end
    end
    %% Plotting the points
    r = norm(W(1,:));
    a = 0.5*pi*rand;
    P = zeros(nrows,2);
    P(1,1) = r*cos(a);
    P(1,2) = r*sin(a);
    for i=2:size(W)
        r = norm(W(i,:));
        if (r ~= 0)
            d = norm(W(i,:)-W(i-1,:));
            a = P(i-1,1);
            b = P(i-1,2);
            syms sol;
            sol = solve((r*cos(sol) - a)^2 + (r*sin(sol) - b)^2 == d^2, sol, 'Real', true);
            sol = double(sol);
            if(numel(sol) && sol(1) >0 && sol(1)<pi/2)
                alpha = sol(1);
            elseif (numel(sol) && sol(2) >0 && sol(2)<pi/2)
                    alpha = sol(2);
            else
                    alpha = rand*range([0 pi/2]);
            end 
        else
            alpha = rand*range([0 pi/2]);
        end
        P(i,1) = r*cos(alpha);
        P(i,2) = r*sin(alpha);
        end
        C = num2cell(P,2);
        y = cell2struct(C,{'cost2d'},2);
end