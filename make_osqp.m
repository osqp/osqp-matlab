function make_osqp(varargin)
% Matlab MEX makefile for OSQP.
%
%    MAKE_OSQP(VARARGIN) is a make file for OSQP solver. It
%    builds OSQP and its components from source.
%
%    WHAT is the last element of VARARGIN and cell array of strings,
%    with the following options:
%
%    {}, '' (empty string) or 'all': build all components and link.
%
%    'osqp_mex': builds the OSQP mex interface and the OSQP library
%
%    Additional commands:
%
%    'clean': Delete all compiled files
%    'purge': Delete all compiled files and copied code generation files


if( nargin == 0 )
    what = {'all'};
    verbose = false;
elseif ( nargin == 1 && ismember('-verbose', varargin) )
    what = {'all'};
    verbose = true;
else
    what = varargin{nargin};
    if(isempty(strfind(what, 'all'))        && ...
        isempty(strfind(what, 'osqp_mex'))  && ...
        isempty(strfind(what, 'clean'))     && ...
        isempty(strfind(what, 'purge')))
            fprintf('No rule to make target "%s", exiting.\n', what);
    end

    verbose = ismember('-verbose', varargin);
end

%% Determine where the various files are all located
% Various parts of the build system
[makefile_path,~,~] = fileparts( which( 'make_osqp.m' ) );
osqp_mex_src_dir = fullfile( makefile_path, 'c_sources' );
osqp_mex_build_dir = fullfile( osqp_mex_src_dir, 'build' );
osqp_cg_src_dir = fullfile( osqp_mex_build_dir, 'codegen_src' );
osqp_cg_dest_dir = fullfile( makefile_path, 'codegen', 'sources' );

% Determine where CMake should look for MATLAB
Matlab_ROOT = strrep( matlabroot, '\', '/' );

%% Try to unlock any pre-existing version of osqp_mex
% this prevents compile errors if a user builds, runs osqp
% and then tries to recompile
if(mislocked('osqp_mex'))
    munlock('osqp_mex');
end

%% Configure, build and install the OSQP mex interface
if( any(strcmpi(what,'osqp_mex')) || any(strcmpi(what,'all')) )
   fprintf('Compiling OSQP solver mex interface...\n');

    % Create build for the mex file and go inside
    if exist( osqp_mex_build_dir, 'dir' )
        rmdir( osqp_mex_build_dir, 's' );
    end
    mkdir( osqp_mex_build_dir );
%    cd( osqp_mex_build_dir );

    % Extend path for CMake mac (via Homebrew)
    PATH = getenv('PATH');
    if( (ismac) && (isempty(strfind(PATH, '/usr/local/bin'))) )
        setenv('PATH', [PATH ':/usr/local/bin']);
    end

    %% Configure CMake for the mex interface
    fprintf('  Configuring...' )
    [status, output] = system( sprintf( 'cmake -B %s -S %s -DMatlab_ROOT_DIR=\"%s\"', osqp_mex_build_dir, osqp_mex_src_dir, Matlab_ROOT ), 'LD_LIBRARY_PATH', '' );
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
    [status, output] = system( sprintf( 'cmake --build %s', osqp_mex_build_dir ), 'LD_LIBRARY_PATH', '' );
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
    [err, errmsg, ~] = copyfile( [osqp_mex_build_dir, filesep, 'osqp_mex.mex*'],  makefile_path );
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

%% Clean and purge
if( any(strcmpi(what,'clean')) || any(strcmpi(what,'purge')) )
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

    fprintf('\t\t[done]\n');

    %% Purge only
    if( any(strcmpi(what,'purge')) )
        fprintf('Cleaning OSQP codegen directories...');

        % Delete codegen files
        if exist(osqp_cg_dest_dir, 'dir')
            rmdir(osqp_cg_dest_dir, 's');
        end

        fprintf('\t\t\t[done]\n');
    end

end

end
