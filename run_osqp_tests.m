import matlab.unittest.TestSuite;

[osqp_path,~,~] = fileparts(which('osqp.m'));
unittest_dir = fullfile(osqp_path, 'unittests');
suiteFolder = TestSuite.fromFolder(unittest_dir);

% Solve individual test file
%suiteFolder = TestSuite.fromFile('unittests/dual_infeasibility_tests.m');

% Run all suite
result = run(suiteFolder);
