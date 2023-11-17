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
    %   setup             - configure solver with problem data
    %   solve             - solve the QP
    %   update            - modify problem vectors
    %   warm_start        - set warm starting variables x and y
    %
    %   default_settings  - create default settings structure
    %   current_settings  - get the current solver settings structure
    %   update_settings   - update the current solver settings structure
    %
    %   get_dimensions    - get the number of variables and constraints
    %   version           - return OSQP version
    %   constant          - return a OSQP internal constant
    %
    %   codegen           - generate embeddable C code for the problem


    properties (SetAccess = private, Hidden = true)
        objectHandle % Handle to underlying C instance
    end
    methods(Static)
        %%
        function out = default_settings()
            % DEFAULT_SETTINGS get the default solver settings structure
            out = osqp_mex('default_settings', 'static');

	        % Convert linsys solver to string
	        out.linsys_solver = linsys_solver_to_string(out.linsys_solver);

        end
        
        %%
        function out = constant(constant_name)
            % CONSTANT Return solver constant
            %   C = CONSTANT(CONSTANT_NAME) return constant called CONSTANT_NAME
            out = osqp_mex('constant', 'static', constant_name);
        end
        
        %%
        function out = version()
            % Return OSQP version
            out = osqp_mex('version', 'static');
        end
        
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
        function update_settings(this,varargin)
            % UPDATE_SETTINGS update the current solver settings structure

            %second input 'false' means that this is *not* a settings
            %initialization, so some parameter/values will be disallowed
            newSettings = validateSettings(this,false,varargin{:});

            %write the solver settings.  C-mex does not check input
            %data or protect against disallowed parameter modifications
            osqp_mex('update_settings', this.objectHandle, newSettings);

        end

        %%
        function [n,m]  = get_dimensions(this)
            % GET_DIMENSIONS get the number of variables and constraints

            [n,m] = osqp_mex('get_dimensions', this.objectHandle);

        end

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
            theSettings = validateSettings(this,true,varargin{6:end});

            [varargout{1:nargout}] = osqp_mex('setup', this.objectHandle, n,m,P,q,A,l,u,theSettings);

        end


        %%

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


            % Decide which function to call
            if (~isempty(x) && isempty(y))
                osqp_mex('warm_start_x', this.objectHandle, x);
                return;
            end

            if (isempty(x) && ~isempty(y))
                osqp_mex('warm_start_y', this.objectHandle, y);
            end

            if (~isempty(x) && ~isempty(y))
                osqp_mex('warm_start', this.objectHandle, x, y);
            end

            if (isempty(x) && isempty(y))
                error('Unrecognized fields');
            end

        end

        %%
        function varargout = solve(this, varargin)
            % SOLVE solve the QP

            nargoutchk(0,1);  %either return nothing (but still solve), or a single output structure
            [out.x, out.y, out.prim_inf_cert, out.dual_inf_cert, out.info] = osqp_mex('solve', this.objectHandle);
            if(nargout)
                varargout{1} = out;
            end
            return;
        end

        %%
        function codegen(this, target_dir, varargin)
            % CODEGEN generate C code for the parametric problem
            %
            %   codegen(target_dir,options)

            % Parse input arguments
            p = inputParser;
            defaultProject = '';
            expectedProject = {'', 'Makefile', 'MinGW Makefiles', 'Unix Makefiles', 'CodeBlocks', 'Xcode'};
            defaultParams = 'vectors';
            expectedParams = {'vectors', 'matrices'};
            defaultMexname = 'emosqp';
            defaultFloat = false;
            defaultLong = true;
            defaultFW = false;

            addRequired(p, 'target_dir', @isstr);
            addParameter(p, 'project_type', defaultProject, ...
                         @(x) ischar(validatestring(x, expectedProject)));
            addParameter(p, 'parameters', defaultParams, ...
                         @(x) ischar(validatestring(x, expectedParams)));
            addParameter(p, 'mexname', defaultMexname, @isstr);
            addParameter(p, 'FLOAT', defaultFloat, @islogical);
            addParameter(p, 'LONG', defaultLong, @islogical);
            addParameter(p, 'force_rewrite', defaultFW, @islogical);

            parse(p, target_dir, varargin{:});

            % Set internal variables
            if strcmp(p.Results.parameters, 'vectors')
                embedded = 1;
            else
                embedded = 2;
            end
            if p.Results.FLOAT
                float_flag = 'ON';
            else
                float_flag = 'OFF';
            end
            if p.Results.LONG
                long_flag = 'ON';
            else
                long_flag = 'OFF';
            end
            if strcmp(p.Results.project_type, 'Makefile')
                if (ispc)
                    project_type = 'MinGW Makefiles';   % Windows
                elseif (ismac || isunix)
                    project_type = 'Unix Makefiles';    % Unix
                end
            else
                project_type = p.Results.project_type;
            end

            % Check whether the specified directory already exists
            if exist(target_dir, 'dir')
                if p.Results.force_rewrite
                    rmdir(target_dir, 's');
                else
                    while(1)
                        prompt = sprintf('Directory "%s" already exists. Do you want to replace it? y/n [y]: ', target_dir);
                        str = input(prompt, 's');

                        if any(strcmpi(str, {'','y'}))
                            rmdir(target_dir, 's');
                            break;
                        elseif strcmpi(str, 'n')
                            return;
                        end
                    end
                end
            end

            % Import OSQP path
            [osqp_path,~,~] = fileparts(which('osqp.m'));

            % Add codegen directory to path
            addpath(fullfile(osqp_path, 'codegen'));

            % Path of osqp module
            cg_dir = fullfile(osqp_path, 'codegen');
            files_to_generate_path = fullfile(cg_dir, 'files_to_generate');

            % Get workspace structure
            work = osqp_mex('get_workspace', this.objectHandle);

            % Make target directory
            fprintf('Creating target directories...\t\t\t\t\t');
            target_configure_dir = fullfile(target_dir, 'configure');
            target_include_dir = fullfile(target_dir, 'include');
            target_src_dir = fullfile(target_dir, 'src');

            if ~exist(target_dir, 'dir')
                mkdir(target_dir);
            end
            if ~exist(target_configure_dir, 'dir')
                mkdir(target_configure_dir);
            end
            if ~exist(target_include_dir, 'dir')
                mkdir(target_include_dir);
            end
            if ~exist(target_src_dir, 'dir')
                mkdir(fullfile(target_src_dir, 'osqp'));
            end
            fprintf('[done]\n');

            % Copy source files to target directory
            fprintf('Copying OSQP source files...\t\t\t\t\t');
            cdir   = fullfile(cg_dir, 'sources', 'src');
            cfiles = dir(fullfile(cdir, '*.c'));
            for i = 1 : length(cfiles)
                if embedded == 1
                    % Do not copy kkt.c if embedded is 1
                    if ~strcmp(cfiles(i).name, 'kkt.c')
                        copyfile(fullfile(cdir, cfiles(i).name), ...
                            fullfile(target_src_dir, 'osqp', cfiles(i).name));    
                    end
                else
                    copyfile(fullfile(cdir, cfiles(i).name), ...
                        fullfile(target_src_dir, 'osqp', cfiles(i).name));
                end
            end
            configure_dir = fullfile(cg_dir, 'sources', 'configure');
            configure_files = dir(fullfile(configure_dir, '*.h.in'));
            for i = 1 : length(configure_files)
                copyfile(fullfile(configure_dir, configure_files(i).name), ...
                    fullfile(target_configure_dir, configure_files(i).name));
            end
            hdir   = fullfile(cg_dir, 'sources', 'include');
            hfiles = dir(fullfile(hdir, '*.h'));
            for i = 1 : length(hfiles)
                if embedded == 1
                    % Do not copy kkt.h if embedded is 1
                    if ~strcmp(hfiles(i).name, 'kkt.h')
                        copyfile(fullfile(hdir, hfiles(i).name), ...
                            fullfile(target_include_dir, hfiles(i).name));  
                    end
                else
                    copyfile(fullfile(hdir, hfiles(i).name), ...
                        fullfile(target_include_dir, hfiles(i).name));
                end
            end

                % Copy cmake files
            copyfile(fullfile(cdir, 'CMakeLists.txt'), ...
                    fullfile(target_src_dir, 'osqp', 'CMakeLists.txt'));
            copyfile(fullfile(hdir, 'CMakeLists.txt'), ...
                    fullfile(target_include_dir, 'CMakeLists.txt'));
            fprintf('[done]\n');

            % Copy example.c
            copyfile(fullfile(files_to_generate_path, 'example.c'), target_src_dir);

            % Render CMakeLists.txt
            fidi = fopen(fullfile(files_to_generate_path, 'CMakeLists.txt'),'r');
            fido = fopen(fullfile(target_dir, 'CMakeLists.txt'),'w');
            while ~feof(fidi)
                l = fgetl(fidi);   % read line
                % Replace EMBEDDED_FLAG in CMakeLists.txt by a numerical value
                newl = strrep(l, 'EMBEDDED_FLAG', num2str(embedded));
                fprintf(fido, '%s\n', newl);
            end
            fclose(fidi);
            fclose(fido);

            % Render workspace.h and workspace.c
            work_hfile = fullfile(target_include_dir, 'workspace.h');
            work_cfile = fullfile(target_src_dir, 'osqp', 'workspace.c');
            fprintf('Generating workspace.h/.c...\t\t\t\t\t\t');
            render_workspace(work, work_hfile, work_cfile, embedded);
            fprintf('[done]\n');

            % Create project
            if ~isempty(project_type)
	    
                % Extend path for CMake mac (via Homebrew)
                PATH = getenv('PATH');
                if ((ismac) && (isempty(strfind(PATH, '/usr/local/bin'))))
                    setenv('PATH', [PATH ':/usr/local/bin']);
                end
    
                fprintf('Creating project...\t\t\t\t\t\t\t\t');
                orig_dir = pwd;
                cd(target_dir);
                mkdir('build')
                cd('build');
                cmd = sprintf('cmake -G "%s" ..', project_type);
                [status, output] = system(cmd);
                if(status)
                    fprintf('\n');
                    fprintf(output);
                    error('Error configuring CMake environment');
                else
                    fprintf('[done]\n');
                end
                cd(orig_dir);
            end

            % Make mex interface to the generated code
            mex_cfile  = fullfile(files_to_generate_path, 'emosqp_mex.c');
            make_emosqp(target_dir, mex_cfile, embedded, float_flag, long_flag);

            % Rename the mex file
            old_mexfile = ['emosqp_mex.', mexext];
            new_mexfile = [p.Results.mexname, '.', mexext];
            movefile(old_mexfile, new_mexfile);

        end

    end
end



function currentSettings = validateSettings(this,isInitialization,varargin)

%don't allow these fields to be changed
unmodifiableFields = {'scaling', 'linsys_solver'};

%get the current settings
if(isInitialization)
    currentSettings = osqp_mex('default_settings', this.objectHandle);
else
    currentSettings = osqp_mex('current_settings', this.objectHandle);
end

%no settings passed -> return defaults
if(isempty(varargin))
    return;
end

%check for structure style input
if(isstruct(varargin{1}))
    newSettings = varargin{1};
    assert(length(varargin) == 1, 'too many input arguments');
else
    newSettings = struct(varargin{:});
end

%get the osqp settings fields
currentFields = fieldnames(currentSettings);

%get the requested fields in the update
newFields = fieldnames(newSettings);

%check for unknown parameters
badFieldsIdx = find(~ismember(newFields,currentFields));
if(~isempty(badFieldsIdx))
    error('Unrecognized solver setting ''%s'' detected',newFields{badFieldsIdx(1)});
end

%convert linsys_solver string to integer
if ismember('linsys_solver',newFields)
   if ~ischar(newSettings.linsys_solver)
       error('Setting linsys_solver is required to be a string.');
   end
   % Convert linsys_solver to number
    newSettings.linsys_solver = string_to_linsys_solver(newSettings.linsys_solver);
end


%check for disallowed fields if this in not an initialization call
if(~isInitialization)
    badFieldsIdx = find(ismember(newFields,unmodifiableFields));
    for i = badFieldsIdx(:)'
        if(~isequal(newSettings.(newFields{i}),currentSettings.(newFields{i})))
            error('Solver setting ''%s'' can only be changed at solver initialization.', newFields{i});
        end
    end
end


%check that everything is a nonnegative scalar (this check is already
%performed in C)
% for i = 1:length(newFields)
%     val = double(newSettings.(newFields{i}));
%     assert(isscalar(val) & isnumeric(val) & val >= 0, ...
%         'Solver setting ''%s'' not specified as nonnegative scalar', newFields{i});
% end

%everything checks out - merge the newSettings into the current ones
for i = 1:length(newFields)
    currentSettings.(newFields{i}) = double(newSettings.(newFields{i}));
end


end

function [linsys_solver_string] = linsys_solver_to_string(linsys_solver)
% Convert linear systme solver integer to stringh
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


