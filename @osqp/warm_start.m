function warm_start(this, varargin)
    % WARM_START warm start primal and/or dual variables
    %
    %   warm_start('x', x, 'y', y)
    %
    %   or warm_start('x', x)
    %   or warm_start('y', y)


    % Get problem dimensions
    [n, m]  = get_dimensions(this);

    % Get data
    allowedFields = {'x','y'};

    if(isempty(varargin))
        return;
    elseif(length(varargin) == 1)
        if(~isstruct(varargin{1}))
            error('Single input should be a structure with new problem data');
        else
            newData = varargin{1};
        end
    else % param / value style assumed
        newData = struct(varargin{:});
    end

    %check for unknown fields
    newFields = fieldnames(newData);
    badFieldsIdx = find(~ismember(newFields,allowedFields));
    if(~isempty(badFieldsIdx))
         error('Unrecognized input field ''%s'' detected',newFields{badFieldsIdx(1)});
    end

    %get all of the terms.  Nonexistent fields will be passed
    %as empty mxArrays
    try x = double(full(newData.x(:))); catch x = []; end
    try y = double(full(newData.y(:))); catch y = []; end

    % Check dimensions
    assert(isempty(x) || length(x) == n, 'input ''x'' is the wrong size');
    assert(isempty(y) || length(y) == m, 'input ''y'' is the wrong size');

    % Only call when there is a vector to update
    if (~isempty(x) || ~isempty(y))
        osqp_mex('warm_start', this.objectHandle, x, y);
    else
        error('Unrecognized fields');
    end
end