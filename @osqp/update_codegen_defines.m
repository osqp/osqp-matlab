function update_codegen_defines(this, varargin)
    % UPDATE_CODEGEN_DEFINES update the current codegen defines

    % Check for structure style input
    if(isstruct(varargin{1}))
        newSettings = varargin{1};
        assert(length(varargin) == 1, 'too many input arguments');
    else
        newSettings = struct(varargin{:});
    end

    % Write the new codegen defiens. The C-function checks for input
    % validity.
    osqp_mex('update_codegen_defines', this.objectHandle, newSettings);
end

