% Demo showing the usage of OSQP from Matlab and the code generation features.
% This problem is the same one that is presented in the osqp_demo.c file.

% The problem data
P = sparse([4., 1.; 1., 2.]);
q = [1; 1];
A = sparse([1., 1; 1, 0; 0, 1]);
l = [1.0; 0.0; 0.0];
u = [1.0; 0.7; 0.7];

% Create the solver
solver = osqp;
solver.setup(P, q, A, l, u, 'verbose', true)

% Solve the problem using the Matlab solver
mat_results = solver.solve()

% Generate the embedded code with the default options
solver.codegen('osqp_demo')

% Solve the problem using the generated code
[x, y, status_val, iter, run_time] = emosqp('solve')
