% Convert linear systme solver integer to string
function [linsys_solver_string] = linsys_solver_to_string(linsys_solver)
    switch linsys_solver
        case osqp.constant('OSQP_UNKNOWN_SOLVER')
            linsys_solver_string = 'unknown solver';
        case osqp.constant('OSQP_DIRECT_SOLVER')
            linsys_solver_string = 'direct solver';
        case osqp.constant('OSQP_INDIRECT_SOLVER')
            linsys_solver_string = 'indirect solver';
        otherwise
            error('Unrecognized linear system solver.');
    end
end
