classdef non_cvx_tests < matlab.unittest.TestCase
    %NON_CVX_TESTS Try to solve a non-convex QP

    properties
        solver
    end

    methods(TestMethodSetup)
        function setup_problem(testCase)
            % Create Problem
            P = sparse([2 5; 5 1]);
            q = [3; 4];
            A = sparse([-1. 0; 0 -1; -1 3; 2  5; 3  4]);
            u = [0; 0; -15.; 100; 80];
            l = -1e20 * ones(length(u), 1);

            % Setup solver
            testCase.solver = osqp;
            testCase.solver.setup(P, q, A, l, u, 'verbose', false);

        end
    end

    methods (Test)
        function test_non_cvx_solve(testCase)
            % Solve with OSQP
            results = testCase.solver.solve();

            % Check if they are close
            testCase.verifyEqual(results.info.status_val, osqp.constant('OSQP_NON_CVX'))
            testCase.verifyTrue(isnan(results.info.obj_val))
        end
        
        function test_nan(testCase)
            testCase.verifyTrue(isnan(osqp.constant('OSQP_NAN')))
        end
    end

end
