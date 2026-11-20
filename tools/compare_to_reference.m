function R = compare_to_reference(rundir, refdir)
%COMPARE_TO_REFERENCE  Compare a run directory with the preserved verified run and report the first divergence.
%   R = compare_to_reference('runs/run_...')                      % reference: reference/baseline_2026-11-20/results
%   R = compare_to_reference(rundir, refdir)
% Checks (exact equality is expected: the code is deterministic on the same MATLAB; floats read back from CSV are
% compared to 1e-12 relative, the precision writetable keeps):
%   1. headline expectations: 200/200 qualified frames for both; IP 415 main + 200 verification = 615; IPC 362;
%      impact frames 24-43: IP 80 + 20 vs IPC 67
%   2. per-frame counts and qualification of both solvers (trajectories.mat, full precision)
%   3. returned-state accuracy and feasibility per frame (r_d, correction, gaps, element dets, IP slack/complementarity)
%   4. phase totals
%   5. the impact-frame replays (impact_sweep_first_limiting_iterates.csv; reference column names mapped)
%   6. the saved replay diagnoses replay_frameNNN.json (limiting vertex, restrictions, centered, IPC same geometry)
% Writes comparison_to_reference.{md,json} into rundir (never into the reference). Returns R.pass and R.first_divergence.
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root, 'src'));
if nargin < 2 || isempty(refdir), refdir = fullfile(root, 'reference', 'baseline_2026-11-20', 'results'); end
rundir = char(rundir); refdir = char(refdir);
N = load(fullfile(rundir, 'trajectories.mat'), 'Rip', 'Ripc');
O = load(fullfile(refdir, 'trajectories.mat'), 'Rip', 'Ripc');
checks = table('Size', [0 5], 'VariableTypes', {'string', 'string', 'string', 'string', 'logical'}, ...
    'VariableNames', {'group', 'quantity', 'expected', 'observed', 'match'});
first = "";

% ---- 1. headline expectations
s = @(F, f) sum([F.(f)]);
imp = 24:43;
H = {
    "IP qualified frames", "200/200", sprintf('%d/200', sum([N.Rip.frames.qualified]))
    "IPC qualified frames", "200/200", sprintf('%d/200', sum([N.Ripc.frames.qualified]))
    "IP main + verification = total", "415 + 200 = 615", sprintf('%d + %d = %d', s(N.Rip.frames, 'main'), s(N.Rip.frames, 'verification'), s(N.Rip.frames, 'total'))
    "IPC main = total", "362 = 362", sprintf('%d = %d', s(N.Ripc.frames, 'main'), s(N.Ripc.frames, 'total'))
    "impact 24-43: IP main + ver vs IPC main", "80 + 20 vs 67", sprintf('%d + %d vs %d', s(N.Rip.frames(imp), 'main'), s(N.Rip.frames(imp), 'verification'), s(N.Ripc.frames(imp), 'main'))};
for q = 1:size(H, 1)
    [checks, first] = add(checks, first, "expected counts", H{q, 1}, H{q, 2}, H{q, 3}, strcmp(H{q, 2}, H{q, 3}));
end

% ---- 2-3. per frame, both solvers
countf = {'main', 'verification', 'discarded', 'fallback', 'total'};
accf = {'rd_final', 'corr_final_tol', 'min_gap_dhat', 'min_det_rel', 'min_alpha_eff', 'comp_dev', 'slack_consistency_dhat', 'unresolved'};
for solver = ["ip", "ipc"]
    a = N.(field(solver)).frames; b = O.(field(solver)).frames;
    [ok, fd, q] = per_frame(a, b, [countf, {'status', 'qualified'}], 'exact');
    [checks, first] = add(checks, first, "per-frame counts and qualification", solver + ": " + strjoin(string([countf, {'status', 'qualified'}]), ', '), ...
        sprintf('%d frames identical', numel(b)), ifelse(ok, sprintf('%d frames identical', numel(a)), sprintf('first divergence: frame %d, %s', fd, q)), ok);
    [ok, fd, q, md] = per_frame(a, b, accf, 'exact');
    [checks, first] = add(checks, first, "returned-state accuracy and feasibility", solver + ": " + strjoin(string(accf), ', '), ...
        'identical', ifelse(ok, 'identical', sprintf('first divergence: frame %d, %s (max abs diff %.3g)', fd, q, md)), ok);
    feas = all([a.min_gap_dhat] > 0) && all([a.min_det_rel] > 0) && all([a.rd_final] <= 5) && all([a.corr_final_tol] <= 1);
    [checks, first] = add(checks, first, "returned-state accuracy and feasibility", solver + ": all gaps > 0, dets > 0, r_d <= 5 F_w, correction <= tol", ...
        'true', string(feas), feas);
end

% ---- 4. phase totals
pt_ref = readtable(fullfile(refdir, 'phase_totals.csv'), 'TextType', 'string');
pt_new = readtable(fullfile(rundir, 'phase_totals.csv'), 'TextType', 'string');
cols = ["main_systems", "verification_systems", "discarded_systems", "fallback_systems", "total_full_size_solves", "failed_frames"];
for q = 1:height(pt_ref)
    r = pt_new(pt_new.phase == pt_ref.phase(q) & pt_new.solver == pt_ref.solver(q), :);
    e = strjoin(string(pt_ref{q, cols}), '/');
    o = "missing"; if ~isempty(r), o = strjoin(string(r{1, cols}), '/'); end
    [checks, first] = add(checks, first, "phase totals", pt_ref.phase(q) + " " + pt_ref.solver(q) + " (main/ver/discarded/fallback/total/failed)", e, o, e == o);
end

% ---- 5. impact replays
sw_ref = readtable(fullfile(refdir, 'impact_sweep_first_limiting_iterates.csv'), 'TextType', 'string');
sw_new = readtable(fullfile(rundir, 'impact_sweep_first_limiting_iterates.csv'), 'TextType', 'string');
map = ["ip_main", "ip_replay_main"; "ip_ver", "ip_replay_ver"; "ip_central_start_main", "ip_central_restart_main";
    "ip_central_start_ver", "ip_central_restart_ver"; "ipc_from_ip_state_main", "ipc_at_ip_state_main";
    "ipc_own_trajectory_main", "ipc_own_trajectory_main"; "rows_carried_at_start", "rows_carried_at_start";
    "rows_central_at_start", "rows_central_at_start"; "k", "k"; "tau_hat", "tau_hat"; "rows", "rows";
    "vertex", "limiting_vertex"; "g_dhat", "g_dhat"; "s_lambda_over_mu", "s_lambda_over_mu";
    "stiffness_ratio_to_central", "stiffness_ratio_to_central"; "motion_over_gap_carried", "motion_over_gap_carried";
    "slack_bound_carried", "slack_bound_carried"; "accepted_alpha_eff", "accepted_alpha_eff";
    "motion_over_gap_centered", "motion_over_gap_centered_same_vertex"; "slack_bound_centered", "slack_bound_centered";
    "centered_limiting_vertex", "centered_limiting_vertex"; "max_closure_centered", "max_closure_centered";
    "motion_over_gap_ipc", "motion_over_gap_ipc_same_vertex"; "ipc_accepted_same_geometry", "ipc_accepted_same_geometry"];
ok = isequal(sw_ref.frame, sw_new.frame); fd = NaN; qd = ""; md = 0;
for fr = sw_ref.frame'
    if ~ok, break; end
    a = sw_ref(sw_ref.frame == fr, :); b = sw_new(sw_new.frame == fr, :);
    for m = 1:size(map, 1)
        x = a.(map(m, 1)); y = b.(map(m, 2));
        d = abs(x - y); rel = d / max(abs(x), 1);
        if ~(isequaln(x, y) || rel <= 1e-12)
            ok = false; fd = fr; qd = map(m, 1); md = d; break;
        end
    end
end
bind_ok = all(sw_ref.binding(~ismissing(sw_ref.binding) & sw_ref.binding ~= "none: full steps") == ...
    sw_new.limiter(~ismissing(sw_ref.binding) & sw_ref.binding ~= "none: full steps"));
[checks, first] = add(checks, first, "impact replays", sprintf('%d frames x %d quantities (replay counts, central restart, IPC at IP state, limiting iterate, centered, IPC same geometry)', height(sw_ref), size(map, 1)), ...
    'identical (1e-12 rel.)', ifelse(ok, 'identical', sprintf('first divergence: frame %d, %s (abs diff %.3g)', fd, qd, md)), ok);
[checks, first] = add(checks, first, "impact replays", "limiter of each first limiting iterate", 'same', string(bind_ok), bind_ok);
ag = readtable(fullfile(rundir, 'impact_replay_agreement.csv'), 'TextType', 'string');
agok = all(ag.ip_identical) && all(ag.ipc_identical(ag.ipc_available));
[checks, first] = add(checks, first, "impact replays", "replays reproduce their own trajectory frames (this run)", ...
    sprintf('%d/%d IP, %d/%d IPC', height(ag), height(ag), sum(ag.ipc_available), sum(ag.ipc_available)), ...
    sprintf('%d/%d IP, %d/%d IPC', sum(ag.ip_identical), height(ag), sum(ag.ipc_identical), sum(ag.ipc_available)), agok);

% ---- 6. saved replay diagnoses
jf = dir(fullfile(refdir, 'replay_frame*.json'));
for q = 1:numel(jf)
    a = jsondecode(fileread(fullfile(refdir, jf(q).name)));
    nf = fullfile(rundir, jf(q).name);
    if ~exist(nf, 'file')
        [checks, first] = add(checks, first, "saved replay diagnoses", string(jf(q).name), 'present', 'missing', false); continue;
    end
    b = jsondecode(fileread(nf));
    L = b.limiting;
    pairs = {
        'IP replay main', a.ip_replay.main, b.ip_replay.main_replay
        'IP replay verification', a.ip_replay.verification, b.ip_replay.verification_replay
        'IPC at IP state main', a.ipc_replay.main, b.ipc_at_ip_state.main
        'IP central restart main', a.ip_central_start_replay.main, b.ip_central_restart.main
        'iterate k', a.iterate.k, L.k
        'accepted fraction', a.iterate.accepted_alpha_eff, L.recorded.alpha_eff
        'limiting vertex', a.vertex.vertex, L.vertex.vertex
        'g / d_hat', a.vertex.g_dhat, L.vertex.g_dhat
        'lambda / F_w', a.vertex.lambda_Fw, L.vertex.lambda_Fw
        's lambda / mu', a.vertex.s_lambda_over_mu, L.vertex.s_lambda_over_mu
        'lambda / s', a.vertex.lambda_over_s_dhat_over_Fw, L.vertex.lambda_over_s_dhat_over_Fw
        'predicted motion / gap', a.vertex.predicted_motion_over_gap, L.vertex.predicted_motion_over_gap
        'slack bound', a.restrictions.slack, L.carried.slack
        'multiplier bound', a.restrictions.multiplier, L.carried.multiplier
        'ground exact crossing', a.restrictions.ground_exact, L.carried.ground_exact
        'merit accepted', a.restrictions.merit_accepted, L.carried.merit_accepted
        'centered slack bound', a.centered.limits.slack, L.centered.limits.slack
        'centered same-vertex motion / gap', a.centered.same_vertex_predicted_motion_over_gap, L.centered.same_vertex_predicted_motion_over_gap
        'IPC same geometry accepted', a.ipc_same_geometry.accepted_alpha, L.ipc_same_geometry.accepted_alpha
        'IPC same geometry motion / gap', a.ipc_same_geometry.same_vertex_predicted_motion_over_gap, L.ipc_same_geometry.same_vertex_predicted_motion_over_gap};
    okj = true; bad = "";
    for p = 1:size(pairs, 1)
        x = pairs{p, 2}; y = pairs{p, 3};
        if ~(isequaln(x, y) || (isnumeric(x) && isnumeric(y) && ~isempty(y) && abs(x - y) <= 1e-12 * max(abs(x), 1)))
            okj = false; bad = sprintf('%s: %g vs %g', pairs{p, 1}, x, y); break;
        end
    end
    [checks, first] = add(checks, first, "saved replay diagnoses", sprintf('%s (%d quantities)', jf(q).name, size(pairs, 1)), ...
        'identical', ifelse(okj, 'identical', bad), okj);
end

R = struct('pass', all(checks.match), 'first_divergence', first, 'rundir', rundir, 'refdir', refdir);
R.checks = checks;
disp(checks);
if R.pass, fprintf('compare_to_reference: all %d checks match\n', height(checks));
else, fprintf('compare_to_reference: FIRST DIVERGENCE: %s\n', first); end
writetable(checks, fullfile(rundir, 'comparison_to_reference.csv'));
fid = fopen(fullfile(rundir, 'comparison_to_reference.md'), 'w');
fprintf(fid, '# Comparison with the reference\n\nrun: `%s`\nreference: `%s`\n\nresult: **%s**\n\n', rundir, refdir, ...
    ifelse(R.pass, 'all checks match', "first divergence: " + first));
fprintf(fid, '| group | quantity | expected | observed | match |\n|---|---|---|---|---|\n');
for q = 1:height(checks)
    fprintf(fid, '| %s | %s | %s | %s | %d |\n', checks.group(q), checks.quantity(q), checks.expected(q), checks.observed(q), checks.match(q));
end
fclose(fid);
end

function [ok, fd, q, md] = per_frame(a, b, fields, ~)
ok = numel(a) == numel(b); fd = NaN; q = ''; md = 0;
if ~ok, fd = min(numel(a), numel(b)) + 1; q = 'number of frames'; return; end
for f = 1:numel(b)
    for k = 1:numel(fields)
        x = b(f).(fields{k}); y = a(f).(fields{k});
        if ~isequaln(x, y)
            ok = false; fd = f; q = fields{k};
            if isnumeric(x), md = max(abs(x - y)); end
            return;
        end
    end
end
end

function [T, first] = add(T, first, g, q, e, o, m)
T = [T; {string(g), string(q), string(e), string(o), logical(m)}];
if ~m && first == "", first = string(g) + " / " + string(q) + ": expected " + string(e) + ", observed " + string(o); end
end

function f = field(solver)
if solver == "ip", f = 'Rip'; else, f = 'Ripc'; end
end

function v = ifelse(c, a, b)
if c, v = a; else, v = b; end
end
