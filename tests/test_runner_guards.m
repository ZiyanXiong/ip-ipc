function T = test_runner_guards()
%TEST_RUNNER_GUARDS  run_experiment refuses an output path inside reference/ and any existing path, before writing.
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root, fullfile(root, 'tools'));
T = table('Size', [0 3], 'VariableTypes', {'string', 'logical', 'string'}, 'VariableNames', {'check', 'passed', 'detail'});
before = manifest_now(fullfile(root, 'reference'));
cases = {
    "inside reference/ (new subfolder)", fullfile(root, 'reference', 'baseline_2026-11-20', 'results', 'new_run'), 'run_experiment:reference'
    "the reference results folder itself", fullfile(root, 'reference', 'baseline_2026-11-20', 'results'), 'run_experiment:reference'
    "relative path into reference/", fullfile('reference', '..', 'reference', 'x'), 'run_experiment:reference'
    "an existing folder", tempdir, 'run_experiment:exists'};
old = pwd; cd(root); c = onCleanup(@() cd(old));
for q = 1:size(cases, 1)
    id = '';
    try
        run_experiment('OutDir', cases{q, 2});
    catch e
        id = e.identifier;
    end
    T = [T; {cases{q, 1}, strcmp(id, cases{q, 3}), string(sprintf('error id: %s', id))}]; %#ok<AGROW>
end
T = [T; {"reference unchanged by the attempts", isequal(before, manifest_now(fullfile(root, 'reference'))), ...
    "file list, sizes and dates under reference/"}];
disp(T);
if ~all(T.passed), error('test_runner_guards: %d of %d checks failed', sum(~T.passed), height(T)); end
fprintf('test_runner_guards: all %d checks passed\n', height(T));
end

function m = manifest_now(d)
f = dir(fullfile(d, '**', '*'));
f = f(~[f.isdir]);
m = sort(arrayfun(@(z) sprintf('%s|%d|%s', fullfile(z.folder, z.name), z.bytes, z.date), f, 'UniformOutput', false));
end
