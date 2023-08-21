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
%    'osqp': builds the OSQP solver using CMake
%
%    'osqp_mex': builds the OSQP mex interface and links it to the OSQP
%    library
%
%    VARARGIN{1:NARGIN-1} specifies the optional flags passed to the compiler
%
%    Additional commands:
%
%    'clean': delete all object files (.o and .obj)
%    'purge' : same as above, and also delete the mex files.


if( nargin == 0 )
    what = {'all'};
    verbose = false;
elseif ( nargin == 1 && ismember('-verbose', varargin) )
    what = {'all'};
    verbose = true;
else
    what = varargin{nargin};
    if(isempty(strfind(what, 'all'))        && ...
        isempty(strfind(what, 'osqp'))      && ...
        isempty(strfind(what, 'osqp_mex'))  && ...
        isempty(strfind(what, 'clean'))     && ...
        isempty(strfind(what, 'purge')))
            fprintf('No rule to make target "%s", exiting.\n', what);
    end

    verbose = ismember('-verbose', varargin);
end

%% Try to unlock any pre-existing version of osqp_mex

% this prevents compile errors if a user builds, runs osqp
% and then tries to recompile
if(mislocked('osqp_mex'))
    munlock('osqp_mex');
end



%% Basic compile commands

% Get make and mex commands
make_cmd = 'cmake --build .';
mex_cmd = sprintf('mex -O -silent');
mex_libs = '';


% Add arguments to cmake and mex compiler
cmake_args = '-DMATLAB=ON';
mexoptflags = '-DMATLAB';

% Add specific generators for windows linux or mac
if (ispc)
    cmake_args = sprintf('%s %s', cmake_args, '-G "MinGW Makefiles"');
else
    cmake_args = sprintf('%s %s', cmake_args, '-G "Unix Makefiles"');
end

% Pass Matlab root to cmake
Matlab_ROOT = strrep(matlabroot, '\', '/');
cmake_args = sprintf('%s %s%s%s', cmake_args, ...
    '-DMatlab_ROOT_DIR="', Matlab_ROOT, '"');

% Add parameters options to mex and cmake
% CTRLC
if (ispc)
   ut = fullfile(matlabroot, 'extern', 'lib', computer('arch'), ...
                 'mingw64', 'libut.lib');
   mex_libs = sprintf('%s "%s"', mex_libs, ut);
else
   mex_libs = sprintf('%s %s', mex_libs, '-lut');
end
% Shared library loading
if (isunix && ~ismac)
   mex_libs = sprintf('%s %s', mex_libs, '-ldl');
end

% Add large arrays support if computer is 64 bit and a pre-2018 version
% Release R2018a corresponds to Matlab version 9.4
if (~isempty(strfind(computer, '64')) && verLessThan('matlab', '9.4'))
    mexoptflags = sprintf('%s %s', mexoptflags, '-largeArrayDims');
end

%Force Matlab to respect old-style usage of mxGetPr in releases after 2018a,
%which use interleaved complex data.   Note that the -R2017b flag is badly
%named since it indicates that non-interleaved complex data model is being used;
%it is not really specific to the release year
if ~verLessThan('matlab', '9.4')
    mexoptflags = sprintf('%s %s', mexoptflags, '-R2017b');
end


% Set optimizer flag
if (~ispc)
    mexoptflags = sprintf('%s %s', mexoptflags, 'COPTIMFLAGS=''-O3''');
end

% Set library extension
lib_ext = '.a';
lib_name = sprintf('libosqpstatic%s', lib_ext);


% Set osqp directory and osqp_build directory
current_dir = pwd;
[makefile_path,~,~] = fileparts(which('make_osqp.m'));
osqp_dir = fullfile(makefile_path, 'osqp_sources');
osqp_build_dir = fullfile(osqp_dir, 'build');
qdldl_dir = fullfile(osqp_dir, 'algebra', '_common', 'lin_sys', 'qdldl');
cg_sources_dir = fullfile('.', 'codegen', 'sources');

% Include directory
inc_dir = [
    fullfile(sprintf(' -I%s', osqp_dir), 'include', 'public'), ...
    fullfile(sprintf(' -I%s', osqp_dir), 'include', 'private'), ...
    fullfile(sprintf(' -I%s', osqp_dir), 'build', 'include', 'public'), ...
    sprintf(' -I%s', qdldl_dir), ...
    fullfile(sprintf(' -I%s', qdldl_dir), 'qdldl_sources', 'include')];


%% OSQP Solver
if( any(strcmpi(what,'osqp')) || any(strcmpi(what,'all')) )
   fprintf('Compiling OSQP solver...');

    % Create build directory and go inside
    if exist(osqp_build_dir, 'dir')
        rmdir(osqp_build_dir, 's');
    end
    mkdir(osqp_build_dir);
    cd(osqp_build_dir);

    % Extend path for CMake mac (via Homebrew)
    PATH = getenv('PATH');
    if ((ismac) && (isempty(strfind(PATH, '/usr/local/bin'))))
        setenv('PATH', [PATH ':/usr/local/bin']);
    end

    % Compile static library with CMake
    [status, output] = system(sprintf('%s %s ..', 'cmake', cmake_args));
    if(status)
        fprintf('\n');
        disp(output);
        error('Error configuring CMake environment');
    elseif(verbose)
        fprintf('\n');
        disp(output);
    end

    [status, output] = system(sprintf('%s %s', make_cmd, '--target osqpstatic'));
    if (status)
        fprintf('\n');
        disp(output);
        error('Error compiling OSQP');
    elseif(verbose)
        fprintf('\n');
        disp(output);
    end


    % Change directory back to matlab interface
    cd(makefile_path);

    % Copy static library to current folder
    lib_origin = fullfile(osqp_build_dir, 'out', lib_name);
    copyfile(lib_origin, lib_name);

    fprintf('\t\t\t\t\t\t[done]\n');

end

%% osqpmex
if( any(strcmpi(what,'osqp_mex')) || any(strcmpi(what,'all')) )
    % Compile interface
    fprintf('Compiling and linking osqpmex...');

    % Compile command
    %cmd = sprintf('%s %s %s %s osqp_mex.cpp', mex_cmd, mexoptflags, inc_dir, lib_name);
    cmd = sprintf('%s %s %s %s osqp_mex.cpp %s', ...
        mex_cmd, mexoptflags, inc_dir, lib_name, mex_libs);

    % Compile
    eval(cmd);
    fprintf('\t\t\t\t\t[done]\n');

end


%% codegen
if( any(strcmpi(what,'codegen')) || any(strcmpi(what,'all')) )
    fprintf('Copying source files for codegen...');

    % Copy C files
    cg_src_dir = fullfile(cg_sources_dir, 'src');
    if ~exist(cg_src_dir, 'dir')
        mkdir(cg_src_dir);
    end
    cdirs  = {fullfile(osqp_dir, 'src'),...
              fullfile(qdldl_dir),...
              fullfile(qdldl_dir, 'qdldl_sources', 'src')};
    for j = 1:length(cdirs)
        cfiles = dir(fullfile(cdirs{j},'*.c'));
        for i = 1 : length(cfiles)
            if ~any(strcmp(cfiles(i).name, {'cs.c', 'ctrlc.c', 'lin_sys.c', 'polish.c'}))
                copyfile(fullfile(cdirs{j}, cfiles(i).name), ...
                    fullfile(cg_src_dir, cfiles(i).name));
            end
        end
    end

    % Copy H files
    cg_include_dir = fullfile(cg_sources_dir, 'include');
    if ~exist(cg_include_dir, 'dir')
        mkdir(cg_include_dir);
    end
    hdirs  = {fullfile(osqp_dir, 'include'),...
              fullfile(qdldl_dir),...
              fullfile(qdldl_dir, 'qdldl_sources', 'include')};
    for j = 1:length(hdirs)
        hfiles = dir(fullfile(hdirs{j},'*.h'));
        for i = 1 : length(hfiles)
            if ~any(strcmp(hfiles(i).name, {'qdldl_types.h', 'osqp_configure.h', ...
                    'cs.h', 'ctrlc.h', 'lin_sys.h', 'polish.h'}))
                copyfile(fullfile(hdirs{j}, hfiles(i).name), ...
                    fullfile(cg_include_dir, hfiles(i).name));
            end
        end
    end

    % Copy configure files
    cg_configure_dir = fullfile(cg_sources_dir, 'configure');
    if ~exist(cg_configure_dir, 'dir')
        mkdir(cg_configure_dir);
    end
    configure_dirs  = {fullfile(osqp_dir, 'configure'),...
                       fullfile(qdldl_dir, 'qdldl_sources', 'configure')};
    for j = 1:length(configure_dirs)
        configure_files = dir(fullfile(configure_dirs{j},'*.h.in'));
        for i = 1 : length(configure_files)
            copyfile(fullfile(configure_dirs{j}, configure_files(i).name), ...
                fullfile(cg_configure_dir, configure_files(i).name));
        end
    end
    
    % Copy cmake files    
    copyfile(fullfile(osqp_dir, 'src', 'CMakeLists.txt'), ...
                    fullfile(cg_src_dir, 'CMakeLists.txt'));
    copyfile(fullfile(osqp_dir, 'include', 'CMakeLists.txt'), ...
                    fullfile(cg_include_dir, 'CMakeLists.txt'));
                
    fprintf('\t\t\t\t\t[done]\n');

end


%% clean
if( any(strcmpi(what,'clean')) || any(strcmpi(what,'purge')) )
    fprintf('Cleaning mex files and library...');

    % Delete mex file
    mexfiles = dir(['*.', mexext]);
    for i = 1 : length(mexfiles)
        delete(mexfiles(i).name);
    end

    % Delete static library
    lib_full_path = fullfile(makefile_path, lib_name);
    if( exist(lib_full_path,'file') )
        delete(lib_full_path);
    end

    fprintf('\t\t\t[done]\n');
end


%% purge
if( any(strcmpi(what,'purge')) )
    fprintf('Cleaning OSQP build and codegen directories...');

    % Delete OSQP build directory
    if exist(osqp_build_dir, 'dir')
        rmdir(osqp_build_dir, 's');
    end

    % Delete codegen files
    if exist(cg_sources_dir, 'dir')
        rmdir(cg_sources_dir, 's');
    end

    fprintf('\t\t[done]\n');
end


%% Go back to the original directory
cd(current_dir);

end
