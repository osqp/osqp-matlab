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