%%
function codegen(this, out, varargin)
    % CODEGEN generate C code for the parametric problem
    %
    %   codegen(target_dir,options)

    % Parse input arguments
    p = inputParser;
    defaultPrefix = 'prob1_';           % Prefix for filenames and C variables; useful if generating multiple problems
    defaultForceRewrite = true;         % Force rewrite if output folder exists?
    defaultParameters = 'vectors';      % What do we wish to update in the generated code?
                                        % One of 'vectors' (allowing update of q/l/u through prob.update_data_vec)
                                        % or 'matrices' (allowing update of P/A/q/l/u
                                        % through prob.update_data_vec or prob.update_data_mat)
    defaultUseFloat = false;            % Use single precision in generated code?
    defaultPrintingEnable = false;      % Enable solver printing?
    defaultProfilingEnable = false;     % Enable solver profiling?
    defaultInterruptEnable = false;     % Enable user interrupt (Ctrl-C)?
    defaultEnableDerivatives = false;   % Enable derivatives?

    addRequired(p, 'out', @isstr);
    addOptional(p, 'prefix', defaultPrefix, @isstr);
    addParameter(p, 'force_rewrite', defaultForceRewrite, @isboolean);
    addParameter(p, 'parameters', defaultParameters, @isstr);
    addParameter(p, 'float_type', defaultUseFloat, @isboolean);
    addParameter(p, 'printing_enable', defaultPrintingEnable, @isboolean);
    addParameter(p, 'profiling_enable', defaultProfilingEnable, @isboolean);
    addParameter(p, 'interrupt_enable', defaultInterruptEnable, @isboolean);
    addParameter(p, 'derivatives_enable', defaultEnableDerivatives, @isboolean);

    parse(p, out, varargin{:});

    % Set internal variables
    if strcmp(p.Results.parameters, 'vectors')
        embedded = 1;
    else
        embedded = 2;
    end


    % Check whether the specified directory already exists
    if exist(out, 'dir')
        while(1)
            prompt = sprintf('Directory "%s" already exists. Do you want to replace it? y/n [y]: ', out);
            str = input(prompt, 's');

            if any(strcmpi(str, {'','y'}))
                rmdir(out, 's');
                break;
            elseif strcmpi(str, 'n')
                return;
            end
        end
    end

    % Import OSQP path
    [osqp_path,~,~] = fileparts(which('osqp.m'));

    % Add codegen directory to path
    addpath(fullfile(osqp_path, 'codegen'));

    % Path of osqp module
    cg_dir = fullfile(osqp_path, '..', 'codegen');
    files_to_generate_path = fullfile(cg_dir, 'files_to_generate');

    % Make target directory
    fprintf('Creating target directories...\t\t\t\t\t');
    target_configure_dir = fullfile(out, 'configure');
    target_include_dir = fullfile(out, 'include');
    target_src_dir = fullfile(out, 'src');

    if ~exist(out, 'dir')
        mkdir(out);
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

    %TODO: Fix the copying stuff
    % % Copy source files to target directory
    % fprintf('Copying OSQP source files...\t\t\t\t\t');
    % cdir   = fullfile(cg_dir, 'sources', 'src');
    % cfiles = dir(fullfile(cdir, '*.c'));
    % for i = 1 : length(cfiles)
    %     if embedded == 1
    %         % Do not copy kkt.c if embedded is 1
    %         if ~strcmp(cfiles(i).name, 'kkt.c')
    %             copyfile(fullfile(cdir, cfiles(i).name), ...
    %                 fullfile(target_src_dir, 'osqp', cfiles(i).name));    
    %         end
    %     else
    %         copyfile(fullfile(cdir, cfiles(i).name), ...
    %             fullfile(target_src_dir, 'osqp', cfiles(i).name));
    %     end
    % end
    % configure_dir = fullfile(cg_dir, 'sources', 'configure');
    % configure_files = dir(fullfile(configure_dir, '*.h.in'));
    % for i = 1 : length(configure_files)
    %     copyfile(fullfile(configure_dir, configure_files(i).name), ...
    %         fullfile(target_configure_dir, configure_files(i).name));
    % end
    % hdir   = fullfile(cg_dir, 'sources', 'inc');
    % hfiles = dir(fullfile(hdir, '*.h'));
    % for i = 1 : length(hfiles)
    %     if embedded == 1
    %         % Do not copy kkt.h if embedded is 1
    %         if ~strcmp(hfiles(i).name, 'kkt.h')
    %             copyfile(fullfile(hdir, hfiles(i).name), ...
    %                 fullfile(target_include_dir, hfiles(i).name));  
    %         end
    %     else
    %         copyfile(fullfile(hdir, hfiles(i).name), ...
    %             fullfile(target_include_dir, hfiles(i).name));
    %     end
    % end
    % 
    %     % Copy cmake files
    % copyfile(fullfile(cdir, 'CMakeLists.txt'), ...
    %         fullfile(target_src_dir, 'osqp', 'CMakeLists.txt'));
    % copyfile(fullfile(hdir, 'CMakeLists.txt'), ...
    %         fullfile(target_include_dir, 'CMakeLists.txt'));
    % fprintf('[done]\n');
    % 
    % % Copy example.c
    % copyfile(fullfile(files_to_generate_path, 'example.c'), target_src_dir);

    % Update codegen defines
    update_codegen_defines(this, 'embedded_mode', embedded, 'float_type', p.Results.float_type, 'printing_enable', p.Results.printing_enable, 'profiling_enable', p.Results.profiling_enable, 'interrupt_enable', p.Results.interrupt_enable, 'derivatives_enable', p.Results.derivatives_enable);
    % Call codegen
    osqp_mex('codegen', this.objectHandle, out, p.Results.prefix);

    % TODO: Do we want to keep this?
    % % Make mex interface to the generated code
    % mex_cfile  = fullfile(files_to_generate_path, 'emosqp_mex.c');
    % make_emosqp(out, mex_cfile, embedded, float_flag, long_flag);
    % 
    % % Rename the mex file
    % old_mexfile = ['emosqp_mex.', mexext];
    % new_mexfile = [p.Results.mexname, '.', mexext];
    % movefile(old_mexfile, new_mexfile);
end