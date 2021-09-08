function make_emosqp(target_dir, mex_cfile, EMBEDDED_FLAG, FLOAT_FLAG, LONG_FLAG)
% Matlab MEX makefile for code generated solver.


% Get make and mex commands
mex_cmd = sprintf('mex -O -silent');

% Add arguments to mex compiler
mexoptflags = '-DMATLAB';

% If running on linux, include the c99 option so GCC uses c99 to compile.
% Otherwise it will throw errors about the comment style
if ( ~ismac() && isunix() )
    mexoptflags = sprintf('%s CFLAGS="$CFLAGS -std=c99"', mexoptflags);
end

% Add embedded flag
cmake_args = sprintf('-DEMBEDDED:INT=%i', EMBEDDED_FLAG);

% Add float flag
cmake_args = sprintf('%s -DDFLOAT:BOOL=%s', cmake_args, FLOAT_FLAG);

% Add long flag
cmake_args = sprintf('%s -DDLONG:BOOL=%s', cmake_args, LONG_FLAG);


% Generate osqp_configure.h file by running cmake
current_dir = pwd;
build_dir = fullfile(target_dir, 'build');
cd(target_dir);
if exist(build_dir, 'dir')
    rmdir('build', 's');
end
mkdir('build');
cd('build');


% Add specific generators for windows linux or mac
if (ispc)
    [status, output] = system(sprintf('%s %s -G "MinGW Makefiles" ..', 'cmake', cmake_args));
else
    [status, output] = system(sprintf('%s %s -G "Unix Makefiles" ..', 'cmake', cmake_args));
end
if(status)
    fprintf('\n');
    disp(output);
    error('Error generating osqp_configure.h');
end
cd(current_dir);

% Set optimizer flag
if (~ispc)
    mexoptflags = sprintf('%s %s', mexoptflags, 'COPTIMFLAGS=''-O3''');
end

% Include directory
inc_dir = fullfile(sprintf(' -I%s', target_dir), 'include');

% Source files
cfiles = '';
src_files = dir(fullfile(target_dir, 'src', 'osqp', '*c'));
for i = 1 : length(src_files)
    cfiles = sprintf('%s %s', cfiles, ...
        fullfile(target_dir, 'src', 'osqp', src_files(i).name));
end

% Compile interface
fprintf('Compiling and linking osqpmex...');

% Compile command
cmd = sprintf('%s %s %s "%s" %s', mex_cmd, mexoptflags, inc_dir, mex_cfile, cfiles);

% Compile
eval(cmd);
fprintf('\t\t\t\t[done]\n');


end
