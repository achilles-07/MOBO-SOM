function [igd] = calculate_igd(true_pareto, approx_pareto)
% This function calculates the Inverted Generational Distance (IGD) between
% the true Pareto front and the approximated Pareto front
%
% Inputs:
%   true_pareto: a matrix containing the true Pareto front
%   approx_pareto: a matrix containing the approximated Pareto front
%
% Outputs:
%   igd: the IGD between the true Pareto front and the approximated Pareto front

n = size(true_pareto, 1);
m = size(true_pareto, 2);
sum_d = 0;

for i = 1:n
    d_i = inf;
    for j = 1:size(approx_pareto, 1)
        d_ij = sqrt(sum((true_pareto(i,:) - approx_pareto(j,:)).^2));
        d_i = min(d_i, d_ij);
    end
    sum_d = sum_d + d_i;
end

igd = sum_d / n;

end
