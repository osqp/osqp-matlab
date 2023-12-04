import matlab.unittest.TestSuite;

[osqp_classpath,~,~] = fileparts( mfilename( 'fullpath' ) );
unittest_dir = fullfile(osqp_classpath, 'unittests');
suiteFolder = TestSuite.fromFolder(unittest_dir);

% Solve individual test file
%suiteFolder = TestSuite.fromFile('unittests/dual_infeasibility_tests.m');

% Run all suite
result = run(suiteFolder);
