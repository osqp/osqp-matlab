function update_settings(this, varargin)
    % UPDATE_SETTINGS update the current solver settings structure

    % Check for structure style input
    if(isstruct(varargin{1}))
        newSettings = varargin{1};
        assert(length(varargin) == 1, 'too many input arguments');
    else
        newSettings = struct(varargin{:});
    end

    % Rho update must be handled special
    if( isfield(newSettings, 'rho') )
        osqp_mex('update_rho', this.objectHandle, newSettings.rho);
        newSettings = rmfield(newSettings, 'rho');
    end

    % Second input 'false' means that this is *not* a settings
    % initialization, so some parameter/values will be disallowed
    newSettings = validate_settings(this, false, varargin{:});

    % Write the solver settings.  C-mex does not check input
    % data or protect against disallowed parameter modifications
    osqp_mex('update_settings', this.objectHandle, newSettings);
end