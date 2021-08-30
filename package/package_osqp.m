function package_osqp()
%   Create OSQP matlab interface package

% Get operative system
if ismac
    platform = 'mac';
elseif isunix
    platform = 'linux';
else
    ispc
    platform = 'windows';
end

% Compile OSQP and get version
fprintf('Compiling OSQP solver\n');
fprintf('---------------------\n');
osqp_dir_matlab = fullfile('..');
cur_dir = pwd;
cd(osqp_dir_matlab);
make_osqp purge;
make_osqp;

% Get OSQP version
s = osqp;
version = s.version;
clear s;
cd(cur_dir)


% Create package
fprintf('Creating Matlab OSQP v%s package\n', version);
fprintf('--------------------------------\n');

% Get package name
package_name = sprintf('osqp-%s-matlab-%s64', version, platform);

% Create package directory
fprintf('Creating package directory %s/...\n', package_name);
if exist(package_name, 'dir')
    rmdir(package_name, 's');
end
mkdir(package_name);
fprintf('[done]\n');

% Copying license
fprintf('Copying license...\n');
copyfile(fullfile(osqp_dir_matlab, 'LICENSE'), ...
	 fullfile(package_name));
fprintf('[done]\n');

% Copying folders
fprintf('Copying folders...\n');
folders_to_copy = {'codegen', 'unittests'};
for i = 1:length(folders_to_copy)
    folder = folders_to_copy{i};
    fprintf('  Copying  %s/%s/...\n', package_name, folder);
    copyfile(fullfile(osqp_dir_matlab, folder), ...
        fullfile(package_name, folder));
end
fprintf('[done]\n');

% Copying files
fprintf('Copying files...\n');
files_to_copy = {sprintf('osqp_mex.%s', mexext),...
    'osqp.m', ...
    'run_osqp_tests.m'};
for i=1:length(files_to_copy)
    file = files_to_copy{i};
    fprintf('  Copying  %s/%s...\n', package_name, file);
    copyfile(fullfile(osqp_dir_matlab, file), ...
        fullfile(package_name, file));
end
fprintf('[done]\n');


fprintf('Compressing files to %s.tar.gz\n', package_name);
tar(sprintf('%s.tar.gz', package_name), package_name);
