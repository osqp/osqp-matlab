classdef non_cvx_tests < matlab.unittest.TestCase
    %NON_CVX_TESTS Try to solve a non-convex QP

    properties
        P
        q
        A
        u
        l
        solver
    end

    methods(TestMethodSetup)
        function setup_problem(testCase)
            % Simple QP
            testCase.P = sparse([2 5; 5 1]);
            testCase.q = [3; 4];
            testCase.A = sparse([-1. 0; 0 -1; -1 3; 2  5; 3  4]);
            testCase.u = [0; 0; -15.; 100; 80];
            testCase.l = -1e20 * ones(length(testCase.u), 1);
            testCase.solver = osqp;

        end
    end

    methods (Test)
        function test_non_cvx_small_sigma(testCase)
            try
                % Setup should fail due to (P + sigma I) having a negative eigenvalue
                testCase.solver.setup(testCase.P, testCase.q, ...
                    testCase.A, testCase.l, testCase.u, ...
                    'verbose', false, 'sigma', 1e-6);
                test_setup = 1;
            catch
                test_setup = 0;
            end

            % Assert test_setup flag
            testCase.verifyTrue(test_setup == 0)
        end

        function test_non_cvx_big_sigma(testCase)
            % Setup workspace with new sigma
            testCase.solver.setup(testCase.P, testCase.q, ...
                    testCase.A, testCase.l, testCase.u, ...
                    'verbose', false, 'sigma', 5);

            % Solve problem
            results = testCase.solver.solve();

            % Assert status and objective value
            testCase.verifyEqual(results.info.status_val, ...
                testCase.solver.constant('OSQP_NON_CVX'))
            testCase.verifyTrue(isnan(results.info.obj_val))
        end
        
        function test_nan(testCase)
            testCase.verifyTrue(isnan(osqp.constant('OSQP_NAN')))
        end
    end

end
