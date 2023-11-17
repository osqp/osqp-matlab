%%
function varargout = setup(this, varargin)
    % SETUP configure solver with problem data
    %
    %   setup(P,q,A,l,u,options)

    nargin = length(varargin);

    %dimension checks on user data. Mex function does not
    %perform any checks on inputs, so check everything here
    assert(nargin >= 5, 'incorrect number of inputs');
    [P,q,A,l,u] = deal(varargin{1:5});

    %
    % Get problem dimensions
    %

    % Get number of variables n
    if (isempty(P))
        if (~isempty(q))
            n = length(q);
        else
            if (~isempty(A))
                n = size(A, 2);
            else
                error('The problem does not have any variables');
            end
        end
    else
        n = size(P, 1);
    end

    % Get number of constraints m
    if (isempty(A))
        m = 0;
    else
        m = size(A, 1);
        assert(size(A, 2) == n, 'Incorrect dimension of A');
    end

    %
    % Create sparse matrices and full vectors if they are empty
    %

    if (isempty(P))
        P = sparse(n, n);
    else
        P = sparse(P);
    end
    if (~istriu(P))
        P = triu(P);
    end
    if (isempty(q))
        q = zeros(n, 1);
    else
        q   = full(q(:));
    end

    % Create proper constraints if they are not passed
    if (isempty(A) && (~isempty(l) || ~isempty(u))) || ...
        (~isempty(A) && (isempty(l) && isempty(u)))
        error('A must be supplied together with at least one bound l or u');
    end

    if (~isempty(A) && isempty(l))
        l = -Inf(m, 1);
    end

    if (~isempty(A) && isempty(u))
        u = Inf(m, 1);
    end

    if (isempty(A))
        A = sparse(m, n);
        l = -Inf(m, 1);
        u = Inf(m, 1);
    else
        l  = full(l(:));
        u  = full(u(:));
        A = sparse(A);
    end


    %
    % Check vector dimensions (not checked from the C solver)
    %

    assert(length(q) == n, 'Incorrect dimension of q');
    assert(length(l) == m, 'Incorrect dimension of l');
    assert(length(u) == m, 'Incorrect dimension of u');

    %
    % Convert infinity values to OSQP_INFINITY
    %
    u = min(u, osqp.constant('OSQP_INFTY'));
    l = max(l, -osqp.constant('OSQP_INFTY'));


    %make a settings structure from the remainder of the arguments.
    %'true' means that this is a settings initialization, so all
    %parameter/values are allowed.  No extra inputs will result
    %in default settings being passed back
    theSettings = validate_settings(this,true,varargin{6:end});

    [varargout{1:nargout}] = osqp_mex('setup', this.objectHandle, n,m,P,q,A,l,u,theSettings);
end