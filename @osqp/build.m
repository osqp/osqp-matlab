function build(varargin)
%MAKE_OSQP(target, args...) A make file for OSQP solver to builds the OSQP mex interface.
%
% Build script to build the OSQP mex interfaces.
%
%  'target' is a string that can take one of the following values:
%    empty: Build osqp_mex interface
%    'all': Build all interfaces
%    'osqp_mex': Builds the OSQP mex interface
%    'osqp_cuda': Builds the OSQP mex interface with the CUDA backend
%    'clean': Delete all compiled files
%    'purge': Delete all compiled files and copied code generation files
%
%  Other arguments:
%    'verbose': true for verbose information during build process

    p = inputParser;

    targets = {'all', 'osqp_mex', 'osqp_cuda', 'clean', 'purge'};

    addOptional( p, 'target', 'osqp_mex', @(s) length( validatestring(s, targets) > 0 ) );
    addParameter( p, 'verbose', false );
    addParameter( p, 'toolchain', 'default', @isstring )

    parse( p, varargin{:} );

    target = p.Results.target;
    verbose = p.Results.verbose;
    cuda_toolchain = p.Results.toolchain;

    %% Determine where the various files are all located
    % Various parts of the build system
    [osqp_classpath,~,~] = fileparts( mfilename( 'fullpath' ) );
    osqp_mex_src_dir = fullfile( osqp_classpath, '..', 'c_sources' );
    osqp_mex_build_dir = fullfile( osqp_mex_src_dir, 'build' );
    osqp_mex_cuda_build_dir = fullfile( osqp_mex_src_dir, 'build_cuda' );
    osqp_cg_src_dir = fullfile( osqp_mex_build_dir, 'codegen_src' );
    osqp_cg_dest_dir = fullfile( osqp_classpath, '..', 'codegen', 'sources' );

    % Determine where CMake should look for MATLAB
    Matlab_ROOT = strrep( matlabroot, '\', '/' );

    %% Try to unlock any pre-existing version of osqp_mex
    % this prevents compile errors if a user builds, runs osqp
    % and then tries to recompile
    if(mislocked('osqp_mex'))
        munlock('osqp_mex');
    end

    %% Configure, build and install the OSQP mex interface
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if( any(strcmpi(target, 'osqp_mex')) || any(strcmpi(target, 'all')) )
       fprintf('Compiling OSQP solver mex interface...\n');

        % Create build for the mex file and go inside
        if exist( osqp_mex_build_dir, 'dir' )
            rmdir( osqp_mex_build_dir, 's' );
        end
        mkdir( osqp_mex_build_dir );

        % Extend path for CMake mac (via Homebrew)
        PATH = getenv('PATH');
        if( (ismac) && (isempty(strfind(PATH, '/usr/local/bin'))) )
            setenv('PATH', [PATH ':/usr/local/bin']);
        end

        %% Configure CMake for the mex interface
        fprintf('  Configuring...' )
        [status, output] = system( sprintf( 'cmake -B %s -S %s -DCMAKE_BUILD_TYPE=RelWithDebInfo -DMatlab_ROOT_DIR=\"%s\"', osqp_mex_build_dir, osqp_mex_src_dir, Matlab_ROOT ), 'LD_LIBRARY_PATH', '' );
        if( status )
            fprintf( '\n' );
            disp( output );
            error( 'Error configuring CMake environment' );
        elseif( verbose )
            fprintf( '\n' );
            disp( output );
        else
            fprintf( '\t\t\t\t\t[done]\n' );
        end

        %% Build the mex interface
        fprintf( '  Building...')
        [status, output] = system( sprintf( 'cmake --build %s --config Release', osqp_mex_build_dir ), 'LD_LIBRARY_PATH', '' );
        if( status )
            fprintf( '\n' );
            disp( output );
            error( 'Error compiling OSQP mex interface' );
        elseif( verbose )
            fprintf( '\n' );
            disp( output );
        else
            fprintf( '\t\t\t\t\t\t[done]\n' );
        end

        %% Install various files
        fprintf( '  Installing...' )

        % Copy mex file to root directory for use
        if( ispc )
            [err, errmsg, ~] = copyfile( [osqp_mex_build_dir, filesep, 'Release', filesep, 'osqp_mex.mex*'], [osqp_classpath, filesep, 'private'] );
        else
            [err, errmsg, ~] = copyfile( [osqp_mex_build_dir, filesep, 'osqp_mex.mex*'],  [osqp_classpath, filesep, 'private'] );
        end
        if( ~err )
            fprintf( '\n' )
            disp( errmsg )
            error( '  Error copying mex file' )
        end

        % Copy the code generation source files
        % Create build for the mex file and go inside
        if exist( osqp_cg_dest_dir, 'dir' )
            rmdir( osqp_cg_dest_dir, 's' );
        end
        mkdir( osqp_cg_dest_dir );

        [err, errmsg, ~] = copyfile( [osqp_cg_src_dir, filesep, '*'], osqp_cg_dest_dir );
        if( ~err )
            fprintf( '\n' )
            disp( errmsg )
            error( '  Error copying code generation source files' )
        end

        fprintf( '\t\t\t\t\t\t[done]\n' );
    end

    %% Configure, build and install the OSQP CUDA mex interface
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if( any(strcmpi(target, 'osqp_cuda')) || any(strcmpi(target, 'all')) )
        fprintf('Compiling OSQP CUDA solver mex interface...\n');

        % Create build for the mex file and go inside
        if exist( osqp_mex_cuda_build_dir, 'dir' )
            rmdir( osqp_mex_cuda_build_dir, 's' );
        end
        mkdir( osqp_mex_cuda_build_dir );

        % Extend path for CMake mac (via Homebrew)
        PATH = getenv('PATH');
        if( (ismac) && (isempty(strfind(PATH, '/usr/local/bin'))) )
            setenv('PATH', [PATH ':/usr/local/bin']);
        end

        %% Figure out the CUDA toolkit and compiler location
        if( any(strcmpi(cuda_toolchain, 'default')) || any(strcmpi(target, 'builtin')) )

            if( exist('mexcuda') )
                cmd = "mexcuda -n -v NVCCFLAGS='$NVCCFLAGS -allow-unsupported-compiler' " + osqp_mex_src_dir + filesep + "configure" + filesep + "config.cu";
                x = evalc( cmd );

                [tok, match] = regexp(x, 'CUDA_ROOT : ([^\n])*', 'tokens', 'match');
                cuda_root = tok{1}{1};

                [tok, match] = regexp(x, 'CUDA_LIBS : ([^\n])*', 'tokens', 'match');
                cuda_libs = tok{1}{1};

                if( verbose )
                    fprintf( '  Found CUDA_ROOT=%s\n', cuda_root );
                    fprintf( '  Found CUDA_LIBS=%s\n', cuda_libs );
                end
            elseif( strcmpi(target, 'builtin') )
                error( 'Unable to locate the builtin MATLAB CUDA toolchain' );
            elseif( verbose )
                fprintf( '  Using CMake CUDA toolchain detection\n' );
            end
        else
            cuda_root = cuda_toolchain;

            if( verbose )
                fprintf( '\n' );
                fprintf( '  Using CUDA_ROOT=%s\n', cuda_root );
            end
        end

        %% Configure CMake for the mex interface
        fprintf('  Configuring...' )
        [status, output] = system( sprintf( 'cmake -B %s -S %s -DCMAKE_BUILD_TYPE=RelWithDebInfo -DMATLAB_CUDA=ON -DMATLAB_CUDA_ROOT=\"%s\" -DMATLAB_CUDA_LIBS=\"%s\" -DMatlab_ROOT_DIR=\"%s\"', ...
                                            osqp_mex_cuda_build_dir, ...
                                            osqp_mex_src_dir, ...
                                            cuda_root, ...
                                            cuda_libs, ...
                                            Matlab_ROOT ), 'LD_LIBRARY_PATH', '' );
        if( status )
            fprintf( '\n' );
            disp( output );
            error( 'Error configuring CMake environment' );
        elseif( verbose )
            fprintf( '\n' );
            disp( output );
        else
            fprintf( '\t\t\t\t\t[done]\n' );
        end

        %% Build the mex interface
        fprintf( '  Building...')
        [status, output] = system( sprintf( 'cmake --build %s --config Release', osqp_mex_cuda_build_dir ), 'LD_LIBRARY_PATH', '' );
        if( status )
            fprintf( '\n' );
            disp( output );
            error( 'Error compiling OSQP mex interface' );
        elseif( verbose )
            fprintf( '\n' );
            disp( output );
        else
            fprintf( '\t\t\t\t\t\t[done]\n' );
        end

        %% Install various files
        fprintf( '  Installing...' )

        % Copy mex file to root directory for use
        if( ispc )
            [err, errmsg, ~] = copyfile( [osqp_mex_cuda_build_dir, filesep, 'Release', filesep, 'osqp_mex.mex*'], [osqp_classpath, filesep, 'private'] );
        else
            [err, errmsg, ~] = copyfile( [osqp_mex_cuda_build_dir, filesep, 'osqp_mex.mex*'],  [osqp_classpath, filesep, 'private'] );
        end
        if( ~err )
            fprintf( '\n' )
            disp( errmsg )
            error( '  Error copying mex file' )
        end

        fprintf( '\t\t\t\t\t\t[done]\n' );
    end

    %% Clean and purge
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if( any(strcmpi(target, 'clean')) || any(strcmpi(target, 'purge')) )
        fprintf('Cleaning OSQP mex files and build directory...');

        % Delete mex file
        mexfiles = dir(['*.', mexext]);
        for i = 1 : length(mexfiles)
            delete(mexfiles(i).name);
        end

        % Delete OSQP build directory
        if exist(osqp_mex_build_dir, 'dir')
            rmdir(osqp_mex_build_dir, 's');
        end
        if exist(osqp_mex_cuda_build_dir, 'dir')
            rmdir(osqp_mex_cuda_build_dir, 's');
        end

        fprintf('\t\t[done]\n');

        %% Purge only
        if( any(strcmpi(target, 'purge')) )
            fprintf('Cleaning OSQP codegen directories...');

            % Delete codegen files
            if exist(osqp_cg_dest_dir, 'dir')
                rmdir(osqp_cg_dest_dir, 's');
            end

            fprintf('\t\t\t[done]\n');
        end
    end
end
