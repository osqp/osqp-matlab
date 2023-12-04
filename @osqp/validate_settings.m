function currentSettings = validate_settings(this, isInitialization, varargin)
    % Don't allow these fields to be changed
    unmodifiableFields = {'scaling', 'linsys_solver'};
    
    % Get the current settings
    if(isInitialization)
        currentSettings = osqp_mex('default_settings', this.objectHandle);
    else
        currentSettings = osqp_mex('current_settings', this.objectHandle);
    end
    
    % No settings passed -> return defaults
    if(isempty(varargin))
        return;
    end
    
    % Check for structure style input
    if(isstruct(varargin{1}))
        newSettings = varargin{1};
        assert(length(varargin) == 1, 'too many input arguments');
    else
        newSettings = struct(varargin{:});
    end
    
    % Get the osqp settings fields
    currentFields = fieldnames(currentSettings);
    
    % Get the requested fields in the update
    newFields = fieldnames(newSettings);
    
    % Check for unknown parameters
    badFieldsIdx = find(~ismember(newFields,currentFields));
    if(~isempty(badFieldsIdx))
        error('Unrecognized solver setting ''%s'' detected',newFields{badFieldsIdx(1)});
    end
    
    % Convert linsys_solver string to integer
    if ismember('linsys_solver',newFields)
       if ~ischar(newSettings.linsys_solver)
           error('Setting linsys_solver is required to be a string.');
       end
       % Convert linsys_solver to number
        newSettings.linsys_solver = string_to_linsys_solver(newSettings.linsys_solver);
    end
    
    
    % Check for disallowed fields if this in not an initialization call
    if(~isInitialization)
        badFieldsIdx = find(ismember(newFields,unmodifiableFields));
        for i = badFieldsIdx(:)'
            if(~isequal(newSettings.(newFields{i}),currentSettings.(newFields{i})))
                error('Solver setting ''%s'' can only be changed at solver initialization.', newFields{i});
            end
        end
    end
    
    
    % Check that everything is a nonnegative scalar (this check is already
    % performed in C)
    % for i = 1:length(newFields)
    %     val = double(newSettings.(newFields{i}));
    %     assert(isscalar(val) & isnumeric(val) & val >= 0, ...
    %         'Solver setting ''%s'' not specified as nonnegative scalar', newFields{i});
    % end
    
    % Everything checks out - merge the newSettings into the current ones
    for i = 1:length(newFields)
        currentSettings.(newFields{i}) = double(newSettings.(newFields{i}));
    end
end
    