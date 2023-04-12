% Generate Solution Space

 for i=1:500
    for j = 1:(M+k-1)
        if(j<M)
            pareto(i,j) = rand;
        else 
            pareto(i,j) = 0;
        end
    end
end


% ----------------------------------------------------------------------------------

% Get Pareto front values

for i=1:500
    true_pareto(i,:) = fobj(pareto(i,:)');
end


% Extract Generated Solution Matrix
for i=1:size(rep)
    approx_pareto(i,:) = rep(i).CostOriginal';
end

% -----------------------------------------------------------------------------------

% Calculate IGD

igd = calculate_igd(true_pareto, approx_pareto);




