function T = test_diagnostics()
%TEST_DIAGNOSTICS  Focused checks of diagnose_iterate.m / diagnose_impact.m on real hard-IP iterates.
%   (a) a slack-limited step                         frame 26, k = 0 (tau_hat = 1, vertex 5)
%   (b) a shortened step with no decreasing slack    frame 67, k = 0 (multiplier positivity at the final level);
%       frame 65, k = 0: slacks decrease but their bound is >= 1, so no slack row is reported either
%   (c) full steps                                    frame 29, k = 0 (17 rows) and frame 24, k = 0 (no rows)
%   (d) every main IP iterate of frames 1-67, re-formed with the shared rules, gives the solver's own limiter and
%       accepted fraction (includes frame 2, whose terminal step the line search rejected at round-off level)
%   (e) the frame-level result of (b) is consumed by the JSON writer and the plots
% Returns a table of checks; errors if any fails.
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root, 'src'));
P = ipipc_params(); S = scene_build(P);
R = run_trajectory(S, P, 'ip', 67);
assert(R.complete, 'the IP trajectory to frame 67 must be complete');
T = table('Size', [0 3], 'VariableTypes', {'string', 'logical', 'string'}, 'VariableNames', {'check', 'passed', 'detail'});

% (a) slack-limited
it = capture(S, P, R, 26, 0); L = diagnose_iterate(S, P, it); V = L.vertex; C = L.carried;
T = add(T, "a: slack-limited limiter", strcmp(L.limiter, 'slack positivity') && C.agrees_with_solver, ...
    sprintf('limiter %s, accepted %.6g (solver %.6g)', L.limiter, C.accepted_alpha_eff, it.alpha_eff));
T = add(T, "a: limiting slack row and contact fields", V.applicable && V.has_row && V.vertex == 5 && ...
    C.slack_vertex == 5 && C.decreasing_slacks > 0 && abs(C.slack - it.alpha_slack) == 0 && ~isempty(V.slack_bound) ...
    && abs(V.slack_bound - C.slack) < 1e-12 && abs(V.s_lambda_over_mu - 0.0561) < 1e-3, ...
    sprintf('vertex %d, slack bound %.6g, s lam/mu %.4g, motion/gap %.4g', V.vertex, C.slack, V.s_lambda_over_mu, V.predicted_motion_over_gap));
T = add(T, "a: centered direction analysed", ~isempty(L.centered.limiter) && L.centered.limits.slack > C.slack, ...
    sprintf('centered: %s, slack bound %.4g', L.centered.limiter, L.centered.limits.slack));

% (b) shortened, no decreasing slack, no slack row
it = capture(S, P, R, 67, 0); L = diagnose_iterate(S, P, it); V = L.vertex; C = L.carried;
T = add(T, "b: no decreasing slack, no slack row", C.decreasing_slacks == 0 && C.slack == 1 && isnan(C.slack_vertex), ...
    sprintf('decreasing slacks %d, slack fraction %g, slack vertex %g', C.decreasing_slacks, C.slack, C.slack_vertex));
T = add(T, "b: actual limiter reported", strcmp(L.limiter, 'multiplier positivity') && C.agrees_with_solver && ...
    it.alpha_eff < 1 && C.accepted_alpha_eff == it.alpha_eff, ...
    sprintf('limiter %s, accepted %.6g (solver %.6g, %s)', L.limiter, C.accepted_alpha_eff, it.alpha_eff, it.limiter));
T = add(T, "b: slack fields empty, multiplier row filled", V.applicable && V.has_row && V.vertex == it.lim_vertex && ...
    isempty(V.slack_bound) && isempty(V.raw_slack_crossing) && abs(V.multiplier_bound - C.multiplier) < 1e-12, ...
    sprintf('vertex %d, multiplier bound %.6g, slack bound field empty', V.vertex, V.multiplier_bound));
it = capture(S, P, R, 65, 0); L = diagnose_iterate(S, P, it); C = L.carried;
T = add(T, "b: slacks decrease but do not restrict: no slack row", C.decreasing_slacks > 0 && C.slack == 1 && ...
    isnan(C.slack_vertex) && strcmp(L.limiter, 'multiplier positivity') && C.agrees_with_solver, ...
    sprintf('decreasing slacks %d, slack fraction %g, slack vertex %g, limiter %s', C.decreasing_slacks, C.slack, C.slack_vertex, L.limiter));

% (c) full steps
it = capture(S, P, R, 29, 0); L = diagnose_iterate(S, P, it); V = L.vertex;
T = add(T, "c: full step with rows", strcmp(L.limiter, 'none') && L.carried.accepted_alpha_eff == 1 && ...
    L.carried.agrees_with_solver && ~V.applicable && isempty(V.vertex) && isempty(V.g_dhat) && L.rows > 0 && ...
    isempty(L.centered.same_vertex_predicted_motion_over_gap) && isempty(L.ipc_same_geometry.same_vertex_Jdx_dhat), ...
    sprintf('%d rows, limiter %s, contact fields empty: %s', L.rows, L.limiter, V.reason));
it = capture(S, P, R, 24, 0); L = diagnose_iterate(S, P, it);
T = add(T, "c: full step without rows", L.rows == 0 && strcmp(L.limiter, 'none') && L.carried.agrees_with_solver && ...
    ~L.vertex.applicable && height(L.rows_table) == 0 && ok_json(L), sprintf('%d rows, limiter %s', L.rows, L.limiter));

% (d) every main iterate of frames 1-67
n = 0; bad = {}; nprec = 0; nlim = struct();
for f = 1:67
    st = R.starts{f};
    [~, ~, its] = ip_solve_frame(S, P, st.x, st.xt, f, st.hist, true);
    for q = find(strcmp({its.kind}, 'main'))
        it = its(q);
        if strcmp(it.limiter, 'discarded'), continue; end
        L = diagnose_iterate(S, P, it);
        n = n + 1;
        nprec = nprec + it.by_precision;
        key = matlab.lang.makeValidName(it.limiter);
        if isfield(nlim, key), nlim.(key) = nlim.(key) + 1; else, nlim.(key) = 1; end
        if ~L.carried.agrees_with_solver
            bad{end+1} = sprintf('f%d k%d: %s/%g vs %s/%g', f, it.k, L.limiter, L.carried.accepted_alpha_eff, it.limiter, it.alpha_eff); %#ok<AGROW>
        end
    end
end
T = add(T, "d: shared rules reproduce every iterate", isempty(bad) && n > 150, ...
    sprintf('%d iterates (%s), %d accepted by the roundoff rule; mismatches: %s', n, limiter_counts(nlim), nprec, strjoin(bad, ', ')));

% (e) the frame-level result is consumed by the writers and the plots
Ripc = run_trajectory(S, P, 'ipc', 67);
D = diagnose_impact(S, P, R, Ripc, 67, 'test (b)');
P67 = P; P67.frames = 67;
ph = trajectory_phases(R, Ripc, P67);
tmp = tempname; mkdir(tmp); c = onCleanup(@() rmdir(tmp, 's'));
files = make_plots(R, Ripc, ph, D, tmp);
J = D; J.limiting = rmfield(J.limiting, 'rows_table'); J = rmfield(J, {'ip_iters', 'ipc_iters'});
T = add(T, "e: writers and plots consume a non-slack result", D.limiting.found && ok_json(J) && numel(files) == 2 && ...
    all(cellfun(@(z) exist(z, 'file') == 2, files)) && D.ip_replay.identical && D.ipc_replay.identical, ...
    sprintf('limiter %s; replay agreement IP %d, IPC %d; %d plots', D.limiting.limiter, D.ip_replay.identical, ...
    D.ipc_replay.identical, numel(files)));

disp(T);
if ~all(T.passed), error('test_diagnostics: %d of %d checks failed', sum(~T.passed), height(T)); end
fprintf('test_diagnostics: all %d checks passed\n', height(T));
end

function it = capture(S, P, R, f, k)
st = R.starts{f};
[~, ~, its] = ip_solve_frame(S, P, st.x, st.xt, f, st.hist, true);
q = find(strcmp({its.kind}, 'main') & [its.k] == k, 1);
it = its(q);
ref = R.iters([R.iters.frame] == f & strcmp({R.iters.kind}, 'main') & [R.iters.k] == k);
assert(it.alpha_eff == ref.alpha_eff && strcmp(it.limiter, ref.limiter), 'capture of frame %d k %d differs from the trajectory', f, k);
end

function T = add(T, name, ok, detail)
T = [T; {name, logical(ok), string(detail)}];
end

function tf = ok_json(s)
try
    jsonencode(rmfield_if(s, 'rows_table'), 'ConvertInfAndNaN', true); tf = true;
catch
    tf = false;
end
end

function s = rmfield_if(s, f)
if isfield(s, f), s = rmfield(s, f); end
end

function t = limiter_counts(s)
f = fieldnames(s);
t = strjoin(cellfun(@(k) sprintf('%s %d', k, s.(k)), f, 'UniformOutput', false), ', ');
end
