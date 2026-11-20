function run_tests()
%RUN_TESTS  The focused checks: diagnostics on real iterates, the failed-frame path, and the reference manifest.
%   matlab -batch "run_tests"
% None of them writes into reference/ or runs/ (the plot check uses a temporary folder).
root = fileparts(mfilename('fullpath'));
addpath(fullfile(root, 'src'), fullfile(root, 'tests'), fullfile(root, 'tools'));
test_diagnostics();
test_failure_handling();
test_runner_guards();
assert(verify_reference(), 'run_tests: the reference does not match its manifest');
fprintf('run_tests: all passed\n');
end
