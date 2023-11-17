%%
function update(this,varargin)
    % UPDATE modify the linear cost term and/or lower and upper bounds

    %second input 'false' means that this is *not* a settings
    %initialization, so some parameter/values will be disallowed
    allowedFields = {'q','l','u','Px','Px_idx','Ax','Ax_idx'};

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
    try q = double(full(newData.q(:))); catch q = []; end
    try l = double(full(newData.l(:))); catch l = []; end
    try u = double(full(newData.u(:))); catch u = []; end
    try Px = double(full(newData.Px(:))); catch Px = []; end
    try Px_idx = double(full(newData.Px_idx(:))); catch Px_idx = []; end
    try Ax = double(full(newData.Ax(:))); catch Ax = []; end
    try Ax_idx = double(full(newData.Ax_idx(:))); catch Ax_idx = []; end

    [n,m]  = get_dimensions(this);

    assert(isempty(q) || length(q) == n, 'input ''q'' is the wrong size');
    assert(isempty(l) || length(l) == m, 'input ''u'' is the wrong size');
    assert(isempty(u) || length(u) == m, 'input ''l'' is the wrong size');
    assert(isempty(Px) || isempty(Px_idx) || length(Px) == length(Px_idx), ...
        'inputs ''Px'' and ''Px_idx'' must be the same size');
    assert(isempty(Ax) || isempty(Ax_idx) || length(Ax) == length(Ax_idx), ...
        'inputs ''Ax'' and ''Ax_idx'' must be the same size');

    % Adjust index of Px_idx and Ax_idx to match 0-based indexing
    % in C
    if (~isempty(Px_idx))
        Px_idx = Px_idx - 1;
    end
    if (~isempty(Ax_idx))
        Ax_idx = Ax_idx - 1;
    end
    
    % Convert infinity values to OSQP_INFTY
    if (~isempty(u))
        u = min(u, osqp.constant('OSQP_INFTY'));
    end
    if (~isempty(l))
        l = max(l, -osqp.constant('OSQP_INFTY'));
    end

    %write the new problem data.  C-mex does not protect
    %against unknown fields, but will handle empty values
    osqp_mex('update', this.objectHandle, ...
    q, l, u, Px, Px_idx, length(Px), Ax, Ax_idx, length(Ax));
end