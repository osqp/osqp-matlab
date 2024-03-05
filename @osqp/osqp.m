classdef osqp < handle
    % osqp interface class for OSQP solver
    % This class provides a complete interface to the C implementation
    % of the OSQP solver.
    %
    % osqp Properties:
    %   objectHandle - pointer to the C structure of OSQP solver
    %
    % osqp Methods:
    %
    %   setup                   - configure solver with problem data
    %   solve                   - solve the QP
    %   update                  - modify problem vectors
    %   warm_start              - set warm starting variables x and y
    %
    %   default_settings        - create default settings structure
    %   current_settings        - get the current solver settings structure
    %   update_settings         - update the current solver settings structure
    %
    %   get_dimensions          - get the number of variables and constraints
    %   version                 - return OSQP version
    %   constant                - return a OSQP internal constant
    %
    %   update_codegen_defines  - update the current codegen defines
    %   codegen                 - generate embeddable C code for the problem


    properties(SetAccess = private, Hidden = true)
        objectHandle % Handle to underlying C instance
    end

    methods(Static)
        output = build(varargin)

        %%
        function out = default_settings()
            % DEFAULT_SETTINGS get the default solver settings structure
            out = osqp_mex('default_settings');

	        % Convert linsys solver to string
	        out.linsys_solver = linsys_solver_to_string(out.linsys_solver);
        end
        
        %%
        function out = constant(constant_name)
            % CONSTANT Return solver constant
            %   C = CONSTANT(CONSTANT_NAME) return constant called CONSTANT_NAME
            out = osqp_mex('constant', constant_name);
        end
        
        %%
        function out = version()
            % Return OSQP version
            out = osqp_mex('version');
        end
    end

    methods(Access = private)
        currentSettings = validate_settings(this, isInitialization, varargin)
    end

    methods
        %% Constructor - Create a new solver instance
        function this = osqp(varargin)
            % Construct OSQP solver class
            this.objectHandle = osqp_mex('new', varargin{:});
        end

        %% Destructor - destroy the solver instance
        function delete(this)
            % Destroy OSQP solver class
            osqp_mex('delete', this.objectHandle);
        end

        %%
        function out = current_settings(this)
            % CURRENT_SETTINGS get the current solver settings structure
            out = osqp_mex('current_settings', this.objectHandle);

            % Convert linsys solver to string
            out.linsys_solver = linsys_solver_to_string(out.linsys_solver);
        end

        %%
        function [n,m]  = get_dimensions(this)
            % GET_DIMENSIONS get the number of variables and constraints

            [n,m] = osqp_mex('get_dimensions', this.objectHandle);
        end
    end
end