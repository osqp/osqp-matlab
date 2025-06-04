% SPDX-License-Identifier: Apache-2.0

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

    % Path to codegen source
    cg_dir = fullfile(osqp_path, '..', 'codegen', 'sources');
    copyfile(cg_dir, out);

    % Update codegen defines
    update_codegen_defines(this, 'embedded_mode', embedded, 'float_type', p.Results.float_type, 'printing_enable', p.Results.printing_enable, 'profiling_enable', p.Results.profiling_enable, 'interrupt_enable', p.Results.interrupt_enable, 'derivatives_enable', p.Results.derivatives_enable);
    % Call codegen
    osqp_mex('codegen', this.objectHandle, out, p.Results.prefix);

end