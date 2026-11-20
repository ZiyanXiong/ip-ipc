function out = run_experiment(varargin)
%RUN_EXPERIMENT  2D IP-vs-IPC box drop: one entry point that reproduces the verified benchmark.
%   From the project folder:   matlab -batch "run_experiment"
%   Optional:                  run_experiment('OutDir', 'C:\some\new\folder')   (the folder must not exist yet)
% Every invocation writes a NEW directory, by default runs/run_<yyyyMMdd_HHmmss>_<8 hex>. An existing path and any
% path inside reference/ are refused, so a run cannot overwrite the preserved reference or mix files with another run.
% Outputs: run_info.json (status, MATLAB version, code digests), settings, per-frame and per-iteration traces,
% phase totals and the performance comparison (only between two complete trajectories), the impact-frame replays
% (agreement with the trajectories, the first limiting IP iterate of each frame), plots, summary.{md,json}, the full
% MATLAB state (trajectories.mat) and run_log.txt. Settings: src/ipipc_params.m (unchanged from the reference run).
ipr = inputParser;
addParameter(ipr, 'OutDir', '', @(s) ischar(s) || isstring(s));
parse(ipr, varargin{:});
root = fileparts(mfilename('fullpath'));
addpath(fullfile(root, 'src'));
out = prepare_outdir(root, char(ipr.Results.OutDir));
info = struct('status', 'running', 'started', char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')), ...
    'matlab', version, 'output_dir', out, 'code_sha256', code_digests(root));
write_json(fullfile(out, 'run_info.json'), info);
diary(fullfile(out, 'run_log.txt'));
cleanup = onCleanup(@() diary('off'));
fprintf('run_experiment  %s  MATLAB %s\noutput: %s\n', info.started, version, out);
t_all = tic;

P = ipipc_params();
S = scene_build(P);
fprintf('scene: %d vertices, %d triangles, F_w = %.6g, kappa = %.6g (corridor [%.4g, %.4g]), tol_disp = %.3g m\n', ...
    S.N, S.nT, S.F_w, S.kappa, S.kappa_corridor, S.tol_disp);

% ---- derivative checks
DC = check_derivatives(S, P);
disp(DC.barrier);
fprintf('IPC normal stiffness positive over (0, d_hat): %d (min %.3g)\n', DC.barrier_kn_positive, DC.barrier_min_kn_over_band);
fprintf('mechanics: gradient rel. err %.2e, Hessian rel. err %.2e, projected min eig %.3g\n', ...
    DC.mech_grad_relerr, DC.mech_hess_relerr, DC.mech_projected_min_eig);

% ---- both trajectories from the scene start
t0 = tic; Ripc = run_trajectory(S, P, 'ipc'); t_ipc = toc(t0);
t0 = tic; Rip  = run_trajectory(S, P, 'ip');  t_ip  = toc(t0);
fprintf('trajectories: IPC %d/%d frames qualified (%.1f s), IP %d/%d (%.1f s)\n', Ripc.last_qualified_frame, P.frames, ...
    t_ipc, Rip.last_qualified_frame, P.frames, t_ip);
ph = trajectory_phases(Rip, Ripc, P);
C = compare_performance(Rip, Ripc, ph);
fprintf('%s\n', C.verdict);

% ---- impact frames: replays and the first limiting iterate (diagnostics, not simulation work)
qual_ip = [Rip.frames([Rip.frames.qualified]).frame];
imp = ph.impact(ismember(ph.impact, qual_ip));
Dall = cell(numel(imp), 1);
for q = 1:numel(imp)
    Dall{q} = diagnose_impact(S, P, Rip, Ripc, imp(q), 'impact sweep');
end
sel = select_frame(Rip, imp, P);
Dsel = []; Dmin = [];
if ~isempty(sel.frame)
    Dsel = Dall{imp == sel.frame}; Dsel.label = 'selected';
    if sel.min_frame ~= sel.frame
        Dmin = Dall{imp == sel.min_frame}; Dmin.label = 'impact frame with the smallest accepted IP fraction';
    end
    fprintf('first contact frame %d; impact %d-%d; selected frame %d (%s)\n', ph.first_contact, ph.impact(1), ...
        ph.impact(end), sel.frame, sel.rule);
end
SW = impact_sweep_table(Dall);
AG = agreement_table(Dall);
fprintf('impact replays: %d/%d IP frames and %d/%d IPC frames reproduce their trajectory exactly\n', ...
    sum(AG.ip_identical), height(AG), sum(AG.ipc_identical), sum(AG.ipc_available));

% ---- outputs
write_json(fullfile(out, 'settings.json'), settings_record(P, S));
write_table(fullfile(out, 'derivative_checks_barrier.csv'), DC.barrier);
write_json(fullfile(out, 'derivative_checks.json'), rmfield(DC, 'barrier'));
write_table(fullfile(out, 'frames_ip.csv'), struct2table(Rip.frames, 'AsArray', true));
write_table(fullfile(out, 'frames_ipc.csv'), struct2table(Ripc.frames, 'AsArray', true));
write_table(fullfile(out, 'iters_ip.csv'), iter_table(Rip.iters));
write_table(fullfile(out, 'iters_ipc.csv'), iter_table(Ripc.iters));
write_table(fullfile(out, 'phase_totals.csv'), C.table);
write_json(fullfile(out, 'performance_comparison.json'), rmfield(C, 'table'));
write_table(fullfile(out, 'impact_sweep_first_limiting_iterates.csv'), SW);
write_table(fullfile(out, 'impact_replay_agreement.csv'), AG);
for R = {Rip, Ripc}
    r = R{1};
    if ~isempty(r.failure)
        fj = rmfield(r.failure, {'iters', 'returned_state', 'frame_start'});
        fj.last_qualified_frame = r.last_qualified_frame; fj.time_end = r.time_end;
        fj.note = 'the state was not advanced: x, v and multiplier history stay those of the last qualified frame';
        write_json(fullfile(out, sprintf('failure_%s.json', r.solver)), fj);
    end
end
for q = 1:numel(Dall)
    d = Dall{q};
    if ~isempty(Dsel) && d.frame == Dsel.frame, d = Dsel; elseif ~isempty(Dmin) && d.frame == Dmin.frame, d = Dmin; end
    if d.limiting.found || (~isempty(Dsel) && d.frame == Dsel.frame), dump_diag(out, d); end
end
plots = make_plots(Rip, Ripc, ph, Dsel, fullfile(out, 'plots'), Dmin);
Sm = summarize(S, P, Rip, Ripc, ph, sel, C, Dall, SW, AG);
write_json(fullfile(out, 'summary.json'), Sm);
write_summary_md(fullfile(out, 'summary.md'), Sm, C, Dsel, Dmin, Rip, Ripc, ...
    cellfun(@(z) strrep(z, [out filesep], ''), plots, 'UniformOutput', false), SW, AG);
save(fullfile(out, 'trajectories.mat'), 'P', 'S', 'Rip', 'Ripc', 'ph', 'sel', 'C', 'Dsel', 'Dmin', 'Dall', 'DC', '-v7.3');
info.status = ternary(C.comparable, 'complete', 'incomplete trajectory (see performance_comparison.json)');
info.finished = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
info.seconds = toc(t_all);
info.ip = C.ip; info.ipc = C.ipc;
write_json(fullfile(out, 'run_info.json'), info);
fprintf('done in %.1f s; results in %s\n', toc(t_all), out);
end

% ======================================================================================== output directory
function out = prepare_outdir(root, req)
if isempty(req)
    stamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
    tag = char(java.util.UUID.randomUUID().toString());
    out = fullfile(root, 'runs', sprintf('run_%s_%s', stamp, tag(1:8)));
else
    out = req;
    if ~is_absolute(out), out = fullfile(pwd, out); end
end
out = char(java.io.File(out).getCanonicalPath());
ref = char(java.io.File(fullfile(root, 'reference')).getCanonicalPath());
if strcmpi(out, ref) || startsWith(lower(out), lower([ref filesep]))
    error('run_experiment:reference', 'refusing to write inside the preserved reference: %s', out);
end
if exist(out, 'file') || exist(out, 'dir')
    error('run_experiment:exists', 'refusing to write into an existing path (runs never share a directory): %s', out);
end
[ok, msg, msgid] = mkdir(out);
if ~ok || strcmp(msgid, 'MATLAB:MKDIR:DirectoryExists')
    error('run_experiment:mkdir', 'cannot create a fresh output directory %s (%s)', out, msg);
end
end

function tf = is_absolute(p)
tf = ~isempty(regexp(p, '^([A-Za-z]:[\\/]|[\\/])', 'once'));
end

function d = code_digests(root)
files = [{fullfile(root, 'run_experiment.m')}; fullfile(root, 'src', {dir(fullfile(root, 'src', '*.m')).name})'];
d = struct('file', {}, 'sha256', {});
for q = 1:numel(files)
    [~, n, e] = fileparts(files{q});
    rel = [n e]; if q > 1, rel = ['src/' rel]; end
    d(end+1) = struct('file', rel, 'sha256', sha256_file(files{q})); %#ok<AGROW>
end
end

% ======================================================================================== selection
function sel = select_frame(Rip, imp, P)
sel = struct('frame', [], 'alpha', [], 'min_frame', [], 'min_alpha', [], 'rule', 'no impact frame was qualified by IP');
if isempty(imp), return; end
a = [Rip.frames(imp).min_alpha_eff];
[amin, qm] = min(a);
sel.min_frame = imp(qm); sel.min_alpha = amin;
q = find(a <= P.short_step_threshold, 1);
if ~isempty(q)
    sel.frame = imp(q); sel.alpha = a(q);
    sel.rule = sprintf('first impact frame with an accepted IP position fraction <= %.2f (it is %.4g)', P.short_step_threshold, a(q));
else
    sel.frame = sel.min_frame; sel.alpha = amin;
    sel.rule = sprintf('no impact frame has an accepted IP fraction <= %.2f: the frame with the smallest (%.4g)', ...
        P.short_step_threshold, amin);
end
end

% ======================================================================================== tables
function T = impact_sweep_table(Dall)
n = numel(Dall);
cols = {'frame', 'ip_replay_main', 'ip_replay_ver', 'ipc_own_trajectory_main', 'ipc_at_ip_state_main', ...
    'ip_central_restart_main', 'ip_central_restart_ver', 'rows_carried_at_start', 'rows_central_at_start', ...
    'limiting_found', 'k', 'tau_hat', 'rows', 'limiter', 'limiting_vertex', 'source', 'g_dhat', 's_lambda_over_mu', ...
    'stiffness_ratio_to_central', 'motion_over_gap_carried', 'slack_bound_carried', 'multiplier_bound_carried', ...
    'accepted_alpha_eff', 'agrees_with_solver', 'centered_limiter', 'centered_limiting_vertex', 'slack_bound_centered', ...
    'multiplier_bound_centered', 'centered_accepted_alpha_eff', 'motion_over_gap_centered_same_vertex', ...
    'max_closure_centered', 'motion_over_gap_ipc_same_vertex', 'ipc_accepted_same_geometry', 'diagnostic_solves'};
rows = cell(0, numel(cols));
for q = 1:n
    d = Dall{q};
    r = {d.frame, d.ip_replay.main_replay, d.ip_replay.verification_replay, d.ipc_at_ip_state.ipc_own_trajectory_main, ...
        d.ipc_at_ip_state.main, d.ip_central_restart.main, d.ip_central_restart.verification, ...
        d.same_state.rows_carried_at_start, d.same_state.rows_central_at_start, d.limiting.found};
    if d.limiting.found
        L = d.limiting; V = L.vertex; Rc = L.carried; Cz = L.centered; I = L.ipc_same_geometry;
        r = [r, {L.k, L.tau_hat, L.rows, string(L.limiter), nn(V.vertex), string(V.source), nn(V.g_dhat), ...
            nn(V.s_lambda_over_mu), nn(V.stiffness_ratio_lambda_s_to_central), nn(V.predicted_motion_over_gap), ...
            Rc.slack, Rc.multiplier, Rc.accepted_alpha_eff, Rc.agrees_with_solver, string(Cz.limiter), ...
            nn(Cz.limiting_vertex), Cz.limits.slack, Cz.limits.multiplier, Cz.limits.accepted_alpha_eff, ...
            nn(Cz.same_vertex_predicted_motion_over_gap), nn(Cz.max_closure_over_gap), ...
            nn(I.same_vertex_predicted_motion_over_gap), I.accepted_alpha, d.diagnostic_solves}]; %#ok<AGROW>
    else
        r = [r, {NaN, NaN, NaN, "none: every accepted step is the full step", NaN, "", NaN, NaN, NaN, NaN, NaN, NaN, ...
            NaN, false, "", NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, d.diagnostic_solves}]; %#ok<AGROW>
    end
    rows(end+1, :) = r; %#ok<AGROW>
end
T = cell2table(rows, 'VariableNames', cols);
end

function T = agreement_table(Dall)
cols = {'frame', 'ip_main_replay', 'ip_main_trajectory', 'ip_ver_replay', 'ip_ver_trajectory', 'ip_status_replay', ...
    'ip_status_trajectory', 'ip_same_iterations', 'ip_max_abs_diff_alpha', 'ip_max_abs_diff_alpha_eff', ...
    'ip_max_abs_diff_rd_inf', 'ip_max_abs_diff_dx_tol', 'ip_max_abs_diff_min_gap', 'ip_same_limiters', 'ip_identical', ...
    'ipc_available', 'ipc_main_replay', 'ipc_main_trajectory', 'ipc_status_replay', 'ipc_status_trajectory', ...
    'ipc_same_iterations', 'ipc_max_abs_diff_alpha', 'ipc_max_abs_diff_rd_inf', 'ipc_max_abs_diff_dx_tol', 'ipc_identical'};
rows = cell(0, numel(cols));
for q = 1:numel(Dall)
    a = Dall{q}.ip_replay; b = Dall{q}.ipc_replay;
    r = {Dall{q}.frame, a.main_replay, a.main_trajectory, a.verification_replay, a.verification_trajectory, ...
        string(a.status_replay), string(a.status_trajectory), a.same_iteration_count, a.max_abs_diff_alpha, ...
        a.max_abs_diff_alpha_eff, a.max_abs_diff_rd_inf, a.max_abs_diff_dx_tol, a.max_abs_diff_min_gap, a.same_limiters, ...
        a.identical, b.available, b.main_replay, b.main_trajectory, string(b.status_replay), string(b.status_trajectory), ...
        b.same_iteration_count, b.max_abs_diff_alpha, b.max_abs_diff_rd_inf, b.max_abs_diff_dx_tol, b.identical};
    rows(end+1, :) = r; %#ok<AGROW>
end
T = cell2table(rows, 'VariableNames', cols);
end

function T = iter_table(its)
its = rmfield(its, 'cap');
T = struct2table(its, 'AsArray', true);
end

function dump_diag(out, D)
tag = sprintf('replay_frame%03d', D.frame);
if D.limiting.found, write_table(fullfile(out, [tag '_rows.csv']), D.limiting.rows_table); end
write_table(fullfile(out, [tag '_iters_ip.csv']), iter_table(D.ip_iters));
write_table(fullfile(out, [tag '_iters_ipc_at_ip_state.csv']), iter_table(D.ipc_iters));
J = rmfield(D, {'ip_iters', 'ipc_iters'});
if J.limiting.found, J.limiting = rmfield(J.limiting, 'rows_table'); end
write_json(fullfile(out, [tag '.json']), J);
end

function write_table(f, T)
writetable(T, f);
end

function write_json(f, s)
fid = fopen(f, 'w');
fwrite(fid, jsonencode(s, 'PrettyPrint', true, 'ConvertInfAndNaN', true), 'char');
fclose(fid);
end

function r = settings_record(P, S)
r.parameters = P;
keep = {'N', 'nT', 'mbar', 'F_w', 'L', 'd_hat', 'tol_disp', 'tau0', 'tau_final', 'cent_final', 's_floor', 'kappa', ...
    'kappa_corridor', 'kappa_dt2', 'scene_diagonal', 'kappa_mean_mass', 'mu', 'lam', 'dt'};
for q = 1:numel(keep), r.resolved.(keep{q}) = S.(keep{q}); end
r.resolved.total_mass_per_m = sum(S.m);
r.resolved.body_weight_in_Fw = sum(S.m) / S.mbar;
r.resolved.bottom_vertices = S.bottom';
r.resolved.first_corner_vertex = S.corner_first;
r.resolved.tau_final_hat = S.tau_final / (S.F_w * S.L);
r.resolved.kappa_dt2_over_Fw_dhat4 = S.kappa_dt2 * S.d_hat^4 / S.F_w;
r.normalization = {
    'The 2D box is a plane-strain slab of unit depth (1 m): masses in kg per m of depth, energies and forces per m.'
    'Solver energy E(x) = 1/2 (x - xt)^T M (x - xt) + dt^2 Psi_el(x) [+ kappa dt^2 sum B(D) for IPC]; gradients are force * dt^2.'
    'F_w = mean vertex mass * |g| * dt^2 is one vertex''s weight in solver units; the body weighs N F_w (N vertices).'
    'Force-balance tolerances are in F_w (per-vertex weight), so the 5 F_w target means the same thing in 2D and 3D.'
    'Length unit L = d_hat; tau has the unit of lambda * s = F_w L; tau_hat = tau / (F_w L); tau_final = tol_disp / L in tau_hat.'
    'kappa follows the native corridor (proportional to the mean vertex mass), so kappa dt^2 / F_w does not depend on the depth.'
    'The displacement tolerance tol_disp = velocity_tol * dt is an absolute length (5e-4 m) for both solvers.'};
end

% ======================================================================================== summary
function Sm = summarize(S, P, Rip, Ripc, ph, sel, C, Dall, SW, AG) %#ok<INUSL>
Sm = struct();
Sm.phases = ph; Sm.selection = sel;
Sm.completion = struct('ip', C.ip, 'ipc', C.ipc, 'comparable', C.comparable, 'verdict', C.verdict);
Sm.counts = C.phases;   % NaN everywhere unless both trajectories are complete
fq = @(R) R.frames([R.frames.qualified]);
fi = fq(Rip); fc = fq(Ripc);
Sm.accuracy_over_qualified_frames = struct( ...
    'ip_max_rd_final_Fw', maxn([fi.rd_final]), 'ipc_max_rd_final_Fw', maxn([fc.rd_final]), ...
    'ip_max_correction_over_tol', maxn([fi.corr_final_tol]), 'ipc_max_correction_over_tol', maxn([fc.corr_final_tol]), ...
    'ip_min_gap_dhat', minn([fi.min_gap_dhat]), 'ipc_min_gap_dhat', minn([fc.min_gap_dhat]), ...
    'ip_min_det_rel', minn([fi.min_det_rel]), 'ipc_min_det_rel', minn([fc.min_det_rel]), ...
    'ip_max_comp_dev', maxn([fi.comp_dev]), 'ip_max_slack_consistency_dhat', maxn([fi.slack_consistency_dhat]), ...
    'ip_max_unresolved', maxn([fi.unresolved]));
Sm.ip_events = struct('early_scans', sum([fi.early_scans]), 'early_rows', sum([fi.early_rows]), ...
    'early_fallbacks', sum([fi.early_fallbacks]), 'outer_clipped', sum([fi.outer_clipped]), 'split', sum([fi.split]), ...
    'dual_limited', sum([fi.dual_limited]), 'valid_limited', sum([fi.valid_limited]), 'merit_backtracks', sum([fi.backtracks]), ...
    'cert_rejects', sum([fi.cert_rejects]), 'terminal_continue', sum([fi.terminal_continue]), ...
    'redirect_would_arm', sum([fi.redirect_would_arm]), 'discarded_frames', [fi([fi.discarded] > 0).frame], ...
    'rows_carried_total', sum([fi.hist_reused]), 'rows_central_at_start_total', sum([fi.hist_fresh]));
Sm.ipc_events = struct('matched_rejects', sum([fc.matched_rejects]), 'ls_limit', sum([fc.ls_limit]), ...
    'cert_rejects', sum([fc.cert_rejects]));
it = Rip.iters(strcmp({Rip.iters.kind}, 'main'));
lims = unique({it.limiter});
Sm.ip_limiter_counts = struct('limiter', lims, 'count', cellfun(@(l) sum(strcmp({it.limiter}, l)), lims, 'UniformOutput', false));
if C.comparable && ~isempty(ph.impact)
    d = [Rip.frames(ph.impact).main] - [Ripc.frames(ph.impact).main];
    Sm.impact_main_difference_per_frame = d;
    short = ph.impact([Rip.frames(ph.impact).min_alpha_eff] <= P.short_step_threshold);
    Sm.short_step_frames = short;
    Sm.short_step_frames_ip_main = sum([Rip.frames(short).main]);
    Sm.short_step_frames_ipc_main = sum([Ripc.frames(short).main]);
end
Sm.impact_sweep = struct('frames', numel(Dall), 'ip_replay_main', sum(SW.ip_replay_main), 'ip_replay_ver', sum(SW.ip_replay_ver), ...
    'ip_central_restart_main', sum(SW.ip_central_restart_main), 'ip_central_restart_ver', sum(SW.ip_central_restart_ver), ...
    'ipc_at_ip_state_main', sum(SW.ipc_at_ip_state_main), 'ipc_own_trajectory_main', sum(SW.ipc_own_trajectory_main), ...
    'frames_with_limiting_iterate', sum(SW.limiting_found), ...
    'diagnostic_solves', sum(SW.diagnostic_solves), 'replay_solves', sum(cellfun(@(z) z.replay_solves, Dall)));
Sm.replay_agreement = struct('frames', height(AG), 'ip_identical', sum(AG.ip_identical), ...
    'ipc_identical', sum(AG.ipc_identical), 'ipc_available', sum(AG.ipc_available));
Sm.wall_s = struct('ip', sum([fi.wall_s]), 'ipc', sum([fc.wall_s]), 'note', 'per-solver wall time; compared only if C.comparable');
end

function write_summary_md(f, Sm, C, Dsel, Dmin, Rip, Ripc, plots, SW, AG)
fid = fopen(f, 'w');
p = @(varargin) fprintf(fid, varargin{:});
p('# Auto-generated summary (run_experiment.m)\n\n');
p('**Completion.** IP: %d/%d frames qualified', Sm.completion.ip.last_qualified_frame, Sm.completion.ip.frames_requested);
if ~Sm.completion.ip.complete, p(' (first failure: frame %d, %s)', Sm.completion.ip.first_failure, Sm.completion.ip.failure_reason); end
p('. IPC: %d/%d', Sm.completion.ipc.last_qualified_frame, Sm.completion.ipc.frames_requested);
if ~Sm.completion.ipc.complete, p(' (first failure: frame %d, %s)', Sm.completion.ipc.first_failure, Sm.completion.ipc.failure_reason); end
p('.\n\n**Comparison.** %s\n\n', C.verdict);
p('Phases: %s. First contact frame %g.\n\n', Sm.phases.rule, Sm.phases.first_contact);
p('## Newton solves by phase (each solver on its own trajectory; qualified frames only)\n\n');
p('| phase | frames | solver | qualified | complete | main | verification | discarded | fallback | total full-size | failed frames | not run |\n|---|---|---|---|---|---|---|---|---|---|---|---|\n');
T = C.table;
for q = 1:height(T)
    p('| %s | %s | %s | %d/%d | %d | %d | %d | %d | %d | %d | %d %s | %d |\n', T.phase{q}, T.frames{q}, T.solver{q}, ...
        T.qualified_frames(q), T.n_frames(q), T.trajectory_complete(q), T.main_systems(q), T.verification_systems(q), ...
        T.discarded_systems(q), T.fallback_systems(q), T.total_full_size_solves(q), T.failed_frames(q), ...
        T.failed_frame_list{q}, T.frames_not_run(q));
end
if C.comparable && ~isempty(Sm.phases.impact)
    p('\n## Impact frames, each solver on its own trajectory\n\n| frame | IP main | IP ver | IP min accepted fraction | IPC main | IPC min accepted fraction |\n|---|---|---|---|---|---|\n');
    for fr = Sm.phases.impact
        a = Rip.frames(fr); b = Ripc.frames(fr);
        p('| %d | %d | %d | %.4g | %d | %.4g |\n', fr, a.main, a.verification, a.min_alpha_eff, b.main, b.min_alpha_eff);
    end
end
A = Sm.accuracy_over_qualified_frames;
p('\n## Accuracy at the returned states (qualified frames)\n\n');
p('- max ||r_d||_inf: IP %.4g F_w, IPC %.4g F_w (target 5)\n', A.ip_max_rd_final_Fw, A.ipc_max_rd_final_Fw);
p('- max returned correction / tol: IP %.4g (verification system), IPC %.4g (last main system)\n', A.ip_max_correction_over_tol, A.ipc_max_correction_over_tol);
p('- min vertex gap: IP %.4g d_hat, IPC %.4g d_hat; min det F (rel. rest): IP %.4g, IPC %.4g\n', A.ip_min_gap_dhat, A.ipc_min_gap_dhat, A.ip_min_det_rel, A.ipc_min_det_rel);
p('- IP: max |s lambda/tau - 1| %.4g, max |g - s| %.4g d_hat, max unresolved contacts %g\n', A.ip_max_comp_dev, A.ip_max_slack_consistency_dhat, A.ip_max_unresolved);
p('\n## Replay agreement (impact_replay_agreement.csv)\n\n');
p('%d of %d IP replays at IP''s own frame-start states and %d of %d IPC replays at IPC''s own frame-start states reproduce the trajectory frame exactly (counts, statuses, every iteration''s step fractions, residuals, corrections, gaps and limiters).\n', ...
    sum(AG.ip_identical), height(AG), sum(AG.ipc_identical), sum(AG.ipc_available));
p('\n## First limiting IP iterate of every impact frame (impact_sweep_first_limiting_iterates.csv)\n\n');
p(['IP replay and IP central restart start at IP''s frame-start state; "IPC at IP state" is a cross replay at that same state; ' ...
    '"IPC own" is IPC on its own trajectory.\n\n']);
p('| frame | IP replay | IP central restart | IPC at IP state | IPC own | k | tau_hat | limiter | vertex | s lam/mu | motion/gap | slack bound | centered: limiter, vertex, slack bound, accepted | IPC same vertex motion/gap |\n|---|---|---|---|---|---|---|---|---|---|---|---|---|---|\n');
for q = 1:height(SW)
    p('| %d | %d+%d | %d+%d | %d | %d | %g | %.3g | %s | %g | %.3g | %.3g | %.3g | %s, %g, %.3g, %.3g | %.3g |\n', SW.frame(q), ...
        SW.ip_replay_main(q), SW.ip_replay_ver(q), SW.ip_central_restart_main(q), SW.ip_central_restart_ver(q), ...
        SW.ipc_at_ip_state_main(q), SW.ipc_own_trajectory_main(q), SW.k(q), SW.tau_hat(q), SW.limiter(q), ...
        SW.limiting_vertex(q), SW.s_lambda_over_mu(q), SW.motion_over_gap_carried(q), SW.slack_bound_carried(q), ...
        SW.centered_limiter(q), SW.centered_limiting_vertex(q), SW.slack_bound_centered(q), ...
        SW.centered_accepted_alpha_eff(q), SW.motion_over_gap_ipc_same_vertex(q));
end
z = Sm.impact_sweep;
p('\nTotals over %d impact frames: IP replay %d + %d; IP central restart %d + %d; IPC at IP state %d; IPC own trajectory %d. Diagnostic solves %d, replay solves %d.\n', ...
    z.frames, z.ip_replay_main, z.ip_replay_ver, z.ip_central_restart_main, z.ip_central_restart_ver, z.ipc_at_ip_state_main, ...
    z.ipc_own_trajectory_main, z.diagnostic_solves, z.replay_solves);
for D = {Dsel, Dmin}
    d = D{1}; if isempty(d) || ~d.limiting.found, continue; end
    L = d.limiting; V = L.vertex; R = L.carried; Cz = L.centered; I = L.ipc_same_geometry; As = L.assembly;
    p('\n## Frame %d (%s): first limiting iterate\n\n', d.frame, d.label);
    p('- k=%d, tau_hat %.3g, %d rows; solver: accepted %.4g, limiter %s; re-formed with the shared rules: %s, accepted %.4g (agrees: %d)\n', ...
        L.k, L.tau_hat, L.rows, L.recorded.alpha_eff, L.recorded.limiter, R.limiter, R.accepted_alpha_eff, R.agrees_with_solver);
    if V.applicable && V.has_row
        p('- vertex %d (%s): g = %.4g d_hat, s = %.4g d_hat, lambda = %.4g F_w, mu = %.4g F_w L, s lambda/mu = %.4g, lambda/s = %.4g F_w/d_hat (mu/s^2 = %.4g), predicted motion %.4g x gap\n', ...
            V.vertex, V.source, V.g_dhat, V.s_dhat, V.lambda_Fw, V.mu_FwL, V.s_lambda_over_mu, V.lambda_over_s_dhat_over_Fw, ...
            V.central_stiffness_mu_over_s2_dhat_over_Fw, V.predicted_motion_over_gap);
    else
        p('- %s\n', V.reason);
    end
    p('- restrictions: slack %.4g (v%s, %d decreasing), multiplier %.4g (v%s, %d decreasing), element %.4g, ground exact crossing %.4g (v%s), merit accepted %.4g after %d backtracks\n', ...
        R.slack, num2str(R.slack_vertex), R.decreasing_slacks, R.multiplier, num2str(R.multiplier_vertex), R.decreasing_multipliers, ...
        R.element, R.ground_exact, num2str(R.ground_vertex), R.merit_accepted, R.merit_backtracks);
    p('- centered lambda = mu/s: %s at v%s, slack %.4g, multiplier %.4g, accepted %.4g; direction change %.3g (rel.)\n', ...
        Cz.limiter, num2str(Cz.limiting_vertex), Cz.limits.slack, Cz.limits.multiplier, Cz.limits.accepted_alpha_eff, Cz.direction_rel_change);
    p('- IPC at the same geometry: %d active, CCD %.4g (v%s), accepted %.4g\n', I.active, I.ccd, num2str(I.ccd_vertex), I.accepted_alpha);
    p('- assembly: s = g %d, stiffness rel err %.3g, rhs (general) rel err %.3g, KKT dx rel err %.3g\n', As.s_equals_g, ...
        As.stiffness_rel_err, As.rhs_general_rel_err, As.kkt_dx_rel_err);
end
p('\nPlots: %s\n', strjoin(plots, ', '));
fclose(fid);
end

function v = nn(v)
if isempty(v), v = NaN; end
end

function m = maxn(v)
if isempty(v), m = NaN; else, m = max(v); end
end

function m = minn(v)
if isempty(v), m = NaN; else, m = min(v); end
end

function v = ternary(c, a, b)
if c, v = a; else, v = b; end
end
