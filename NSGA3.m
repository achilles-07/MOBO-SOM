%% Define the problem
n = 20; % number of variables
m = 10; % number of objectives
k = 100; % population size
generations = 200; % number of generations

problem = struct();
problem.name = 'DTLZ5';
problem.params = struct();
problem.params.num_vars = n;
problem.params.num_objs = m;
problem.params.var_lb = zeros(1, n);
problem.params.var_ub = ones(1, n);
problem.params.objfun = @(x) dtlz5(x, m);

%% Run NSGA-II
options = nsgaopt('popsize', k, 'generations', generations, 'crossover', 0.9, 'mutation', 1/n, 'plotinterval', 1);
result = nsga2(problem, options);

%% Plot the results
plot(result.pop(:,m+1), result.pop(:,m+2), 'bo');
xlabel('f_1');
ylabel('f_2');
title('Approximated Pareto front for DTLZ5 with 10 objectives');
