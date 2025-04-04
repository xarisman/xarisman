% Define the objective function
BW = 0.05*pi;
c1 = 0.16;
c2 = 0.83;


objective_function = @(p) ga_function(p, BW, c1, c2);

% Define the number of variables (7 in this case)
num_vars = 3;

% Set bounds for each variable in p
% Adjust these bounds based on your requirements
lb = [0.01, 0.01, 0.01];  % Lower bounds for p variables
ub = [1, 1, 1];         % Upper bounds for p variables

% Set up options for the genetic algorithm
options = optimoptions('ga', ...
    'PopulationSize', 100, ...       % Size of the population
    'MaxGenerations', 200, ...      % Maximum number of generations
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

T_triplos(best_p)
