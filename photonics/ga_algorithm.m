% Define the objective function
objective_function = @(p) func_min(p);

% Define the number of variables (7 in this case)
num_vars = 9;

% Set bounds for each variable in p
% Adjust these bounds based on your requirements
lb = [0.05, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.05];  % Lower bounds for p variables
ub = [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5];         % Upper bounds for p variables

% Set up options for the genetic algorithm
options = optimoptions('ga', ...
    'PopulationSize', 50, ...       % Size of the population
    'MaxGenerations', 100, ...      % Maximum number of generations
    'Display', 'iter', ...          % Display iteration information
    'PlotFcn', {@gaplotbestf});     % Plot best objective function value

% Run the genetic algorithm
% 'ga' will attempt to minimize the objective function
[best_p, best_fval] = ga(objective_function, num_vars, [], [], [], [], lb, ub, [], options);

% Display the results
disp('Optimal parameter values (p):');
disp(best_p);
disp('Optimal objective function value:');
disp(best_fval);