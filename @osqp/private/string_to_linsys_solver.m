function [linsys_solver] = string_to_linsys_solver(linsys_solver_string)
    linsys_solver_string = lower(linsys_solver_string);
    switch linsys_solver_string
        case 'unknown solver'
            linsys_solver = osqp.constant('OSQP_UNKNOWN_SOLVER');
        case 'direct solver'
            linsys_solver = osqp.constant('OSQP_DIRECT_SOLVER');
        case 'indirect solver'
            linsys_solver = osqp.constant('OSQP_INDIRECT_SOLVER');
        % Default solver: QDLDL
        case ''
            linsys_solver = osqp.constant('OSQP_DIRECT_SOLVER');
        otherwise
            warning('Linear system solver not recognized. Using default solver OSQP_DIRECT_SOLVER.')
            linsys_solver = osqp.constant('OSQP_DIRECT_SOLVER');
    end
end
    