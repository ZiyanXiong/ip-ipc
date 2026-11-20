function run_experiment()
%RUN_EXPERIMENT  2D IP-vs-IPC box drop: one entry point that reproduces everything.
%   matlab -batch "cd('E:\ZiyanXiong\github\ipipc-matlab'); run_experiment"
% Writes results/ (settings, per-frame and per-iteration traces, phase totals, same-state replay diagnosis,
% plots, summary.md). Reference (read-only): E:\ZiyanXiong\github\libuipc-ip.
root = fileparts(mfilename('fullpath'));
addpath(fullfile(root, 'src'));
out = fullfile(root, 'results');
if ~exist(out, 'dir'), mkdir(out); end
if exist(fullfile(out, 'run_log.txt'), 'file'), delete(fullfile(out, 'run_log.txt')); end
diary(fullfile(out, 'run_log.txt'));
cleanup = onCleanup(@() diary('off'));
fprintf('run_experiment  %s  MATLAB %s\n', char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')), version);
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
fprintf('trajectories: IPC %.1f s, IP %.1f s\n', t_ipc, t_ip);

% ---- phases, selection, same-state diagnosis
ph = phases(Rip, Ripc, P);
sel = select_frame(Rip, ph, P);
fprintf('first contact frame %d; impact %d-%d; selected frame %d (%s)\n', ph.first_contact, ph.impact(1), ph.impact(end), ...
    sel.frame, sel.rule);
% every impact frame is replayed from IP's frame-start state (diagnostics, not simulation work); the selected frame
% and the frame with the smallest accepted fraction are written out in full
Dall = cell(numel(ph.impact), 1);
for q = 1:numel(ph.impact)
    Dall{q} = diagnose_impact(S, P, Rip, Ripc, ph.impact(q), 'impact sweep');
end
Dsel = Dall{ph.impact == sel.frame}; Dsel.label = 'selected';
Dmin = [];
if sel.min_frame ~= sel.frame
    Dmin = Dall{ph.impact == sel.min_frame}; Dmin.label = 'impact frame with the smallest accepted IP fraction';
end
SW = impact_sweep_table(Dall);
fprintf('impact sweep: %d frames with a limiting iterate\n', sum(~isnan(SW.k)));

% ---- outputs
write_json(fullfile(out, 'settings.json'), settings_record(P, S));
write_table(fullfile(out, 'derivative_checks_barrier.csv'), DC.barrier);
write_json(fullfile(out, 'derivative_checks.json'), rmfield(DC, 'barrier'));
write_table(fullfile(out, 'frames_ip.csv'), struct2table(Rip.frames));
write_table(fullfile(out, 'frames_ipc.csv'), struct2table(Ripc.frames));
write_table(fullfile(out, 'iters_ip.csv'), iter_table(Rip.iters));
write_table(fullfile(out, 'iters_ipc.csv'), iter_table(Ripc.iters));
PT = phase_table(Rip, Ripc, ph);
write_table(fullfile(out, 'phase_totals.csv'), PT);
write_table(fullfile(out, 'impact_sweep_first_limiting_iterates.csv'), SW);
for q = 1:numel(Dall)
    d = Dall{q};
    if d.frame == Dsel.frame, d = Dsel; elseif ~isempty(Dmin) && d.frame == Dmin.frame, d = Dmin; end
    if isfield(d, 'vertex') || d.frame == Dsel.frame, dump_diag(out, d); end
end
plots = make_plots(Rip, Ripc, ph, Dsel, fullfile(out, 'plots'), Dmin);
Sm = summarize(S, P, Rip, Ripc, ph, sel, PT, Dsel, Dmin, DC);
Sm.impact_sweep = sweep_totals(SW, Dall);
write_json(fullfile(out, 'summary.json'), Sm);
write_summary_md(fullfile(out, 'summary.md'), Sm, PT, Dsel, Dmin, Rip, Ripc, plots, SW);
starts_ip = Rip.starts; %#ok<NASGU>
save(fullfile(out, 'trajectories.mat'), 'P', 'S', 'Rip', 'Ripc', 'ph', 'sel', 'Dsel', 'Dmin', 'DC', '-v7.3');
fprintf('done in %.1f s; results in %s\n', toc(t_all), out);
end

% ======================================================================================== phases and selection
function ph = phases(Rip, Ripc, P)
c1 = find([Rip.frames.contact], 1); c2 = find([Ripc.frames.contact], 1);
fc = min([c1, c2]);
nf = numel(Rip.frames);
ph.first_contact = fc; ph.first_contact_ip = c1; ph.first_contact_ipc = c2;
ph.free = 1:fc-1;
ph.impact = fc:min(nf, fc + P.impact_frames - 1);
ph.settling = ph.impact(end)+1:nf;
ph.all = 1:nf;
ph.rule = sprintf(['free flight = frames before the first frame in which either solver has a contact pair; impact = that ' ...
    'frame and the next %d; settling = the rest'], P.impact_frames - 1);
end

function sel = select_frame(Rip, ph, P)
a = [Rip.frames(ph.impact).min_alpha_eff];
q = find(a <= P.short_step_threshold, 1);
[amin, qm] = min(a);
sel.min_frame = ph.impact(qm); sel.min_alpha = amin;
if ~isempty(q)
    sel.frame = ph.impact(q); sel.alpha = a(q);
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
c = struct('frame', [], 'ip_main', [], 'ip_ver', [], 'ip_central_start_main', [], 'ip_central_start_ver', [], ...
    'ipc_from_ip_state_main', [], 'ipc_own_trajectory_main', [], 'rows_carried_at_start', [], 'rows_central_at_start', [], ...
    'k', [], 'tau_hat', [], 'rows', [], 'vertex', [], 'source', strings(0), 'g_dhat', [], 's_lambda_over_mu', [], ...
    'stiffness_ratio_to_central', [], 'motion_over_gap_carried', [], 'slack_bound_carried', [], 'accepted_alpha_eff', [], ...
    'binding', strings(0), 'motion_over_gap_centered', [], 'slack_bound_centered', [], 'centered_limiting_vertex', [], ...
    'max_closure_centered', [], 'motion_over_gap_ipc', [], 'ipc_accepted_same_geometry', [], 'diagnostic_solves', []);
for q = 1:n
    d = Dall{q};
    c.frame(q, 1) = d.frame;
    c.ip_main(q, 1) = d.ip_replay.main; c.ip_ver(q, 1) = d.ip_replay.verification;
    c.ip_central_start_main(q, 1) = d.ip_central_start_replay.main; c.ip_central_start_ver(q, 1) = d.ip_central_start_replay.verification;
    c.ipc_from_ip_state_main(q, 1) = d.ipc_replay.main; c.ipc_own_trajectory_main(q, 1) = d.ipc_replay.trajectory_main_on_own_trajectory;
    c.rows_carried_at_start(q, 1) = d.same_state.rows_carried_at_start; c.rows_central_at_start(q, 1) = d.same_state.rows_central_at_start;
    c.diagnostic_solves(q, 1) = d.diagnostic_solves;
    if isfield(d, 'vertex')
        V = d.vertex; R = d.restrictions; C = d.centered; I = d.ipc_same_geometry;
        c.k(q, 1) = d.iterate.k; c.tau_hat(q, 1) = d.iterate.tau_hat; c.rows(q, 1) = d.iterate.rows; c.vertex(q, 1) = V.vertex;
        c.source(q, 1) = string(V.source); c.g_dhat(q, 1) = V.g_dhat; c.s_lambda_over_mu(q, 1) = V.s_lambda_over_mu;
        c.stiffness_ratio_to_central(q, 1) = V.stiffness_ratio_lambda_s_to_central;
        c.motion_over_gap_carried(q, 1) = V.predicted_motion_over_gap; c.slack_bound_carried(q, 1) = R.slack;
        c.accepted_alpha_eff(q, 1) = d.iterate.accepted_alpha_eff; c.binding(q, 1) = string(R.binding);
        c.motion_over_gap_centered(q, 1) = C.same_vertex_predicted_motion_over_gap; c.slack_bound_centered(q, 1) = C.limits.slack;
        c.centered_limiting_vertex(q, 1) = C.limiting_vertex; c.max_closure_centered(q, 1) = C.max_closure_over_gap;
        c.motion_over_gap_ipc(q, 1) = I.same_vertex_predicted_motion_over_gap; c.ipc_accepted_same_geometry(q, 1) = I.accepted_alpha;
    else
        c.k(q, 1) = NaN; c.tau_hat(q, 1) = NaN; c.rows(q, 1) = NaN; c.vertex(q, 1) = NaN; c.source(q, 1) = "";
        c.g_dhat(q, 1) = NaN; c.s_lambda_over_mu(q, 1) = NaN; c.stiffness_ratio_to_central(q, 1) = NaN;
        c.motion_over_gap_carried(q, 1) = NaN; c.slack_bound_carried(q, 1) = NaN; c.accepted_alpha_eff(q, 1) = NaN;
        c.binding(q, 1) = "none: full steps"; c.motion_over_gap_centered(q, 1) = NaN; c.slack_bound_centered(q, 1) = NaN;
        c.centered_limiting_vertex(q, 1) = NaN; c.max_closure_centered(q, 1) = NaN; c.motion_over_gap_ipc(q, 1) = NaN;
        c.ipc_accepted_same_geometry(q, 1) = NaN;
    end
end
T = struct2table(c);
end

function s = sweep_totals(SW, Dall)
s.frames = numel(Dall);
s.ip_main = sum(SW.ip_main); s.ip_ver = sum(SW.ip_ver);
s.ip_central_start_main = sum(SW.ip_central_start_main); s.ip_central_start_ver = sum(SW.ip_central_start_ver);
s.ipc_from_ip_state_main = sum(SW.ipc_from_ip_state_main); s.ipc_own_trajectory_main = sum(SW.ipc_own_trajectory_main);
lim = ~isnan(SW.k);
s.frames_with_limiting_iterate = sum(lim);
s.limiting_at_tau1 = sum(lim & SW.tau_hat > 0.5);
s.limiting_slack_binding = sum(lim & SW.binding == "slack positivity");
s.limiting_vertex_carried = sum(lim & contains(SW.source, "carried"));
s.centered_still_closes_beyond_gap = sum(lim & SW.motion_over_gap_centered < -1);
s.centered_slack_bound_ratio = SW.slack_bound_centered(lim) ./ SW.slack_bound_carried(lim);
s.diagnostic_solves = sum(SW.diagnostic_solves);
s.replay_solves = sum(cellfun(@(d) d.replay_solves, Dall));
end

function T = iter_table(its)
its = rmfield(its, 'cap');
T = struct2table(its, 'AsArray', true);
end

function T = phase_table(Rip, Ripc, ph)
names = {'free flight', 'impact', 'settling', 'total'};
sets = {ph.free, ph.impact, ph.settling, ph.all};
rows = {};
for q = 1:4
    fs = sets{q};
    for R = {Rip, Ripc}
        r = R{1}; F = r.frames(fs);
        rows(end+1, :) = {names{q}, sprintf('%d-%d', fs(1), fs(end)), r.solver, numel(fs), sum([F.main]), ...
            sum([F.verification]), sum([F.discarded]), sum([F.fallback]), sum([F.main]) + sum([F.verification]), ...
            sum(~[F.qualified]), strjoin(arrayfun(@(z) sprintf('%d', z), fs(~[F.qualified]), 'UniformOutput', false), ' '), ...
            sum([F.wall_s])}; %#ok<AGROW>
    end
end
T = cell2table(rows, 'VariableNames', {'phase', 'frames', 'solver', 'n_frames', 'main_systems', 'verification_systems', ...
    'discarded_systems', 'fallback_systems', 'total_full_size_solves', 'failed_frames', 'failed_frame_list', 'wall_s'});
end

function dump_diag(out, D)
tag = sprintf('replay_frame%03d', D.frame);
if isfield(D, 'rows_table'), write_table(fullfile(out, [tag '_rows.csv']), D.rows_table); end
write_table(fullfile(out, [tag '_iters_ip.csv']), iter_table(D.ip_iters));
write_table(fullfile(out, [tag '_iters_ipc.csv']), iter_table(D.ipc_iters));
J = rmfield(D, intersect(fieldnames(D), {'rows_table', 'ip_iters', 'ipc_iters'}));
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
function Sm = summarize(S, P, Rip, Ripc, ph, sel, PT, Dsel, Dmin, DC) %#ok<INUSD>
Sm = struct();
Sm.phases = ph; Sm.selection = sel;
get = @(phase, solver, col) PT{strcmp(PT.phase, phase) & strcmp(PT.solver, solver), col};
for q = {'free flight', 'impact', 'settling', 'total'}
    k = matlab.lang.makeValidName(q{1});
    Sm.counts.(k) = struct('ip_main', get(q{1}, 'ip', 'main_systems'), 'ip_verification', get(q{1}, 'ip', 'verification_systems'), ...
        'ip_total', get(q{1}, 'ip', 'total_full_size_solves'), 'ipc_main', get(q{1}, 'ipc', 'main_systems'), ...
        'ipc_total', get(q{1}, 'ipc', 'total_full_size_solves'));
end
% short-step attribution on IP's own trajectory
it = Rip.iters(strcmp({Rip.iters.kind}, 'main'));
inimp = ismember([it.frame], ph.impact);
slack = contains({it.limiter}, 'slack');
Sm.ip_impact_slack_limited_systems = sum(inimp & slack & [it.alpha_eff] < 1);
Sm.ip_impact_short_systems = sum(inimp & [it.alpha_eff] <= P.short_step_threshold & [it.alpha_eff] > 0);
Sm.ip_impact_slack_limited_first_iterate_frames = sum(inimp & slack & [it.k] == 0 & [it.alpha_eff] < 1);
first_it = it(inimp & [it.k] == 0 & [it.rows] > 0);
Sm.ip_impact_first_iterates = numel(first_it);
Sm.ip_impact_first_iterate_alpha = [first_it.alpha_eff];
Sm.ip_impact_first_iterate_slam_median = [first_it.slam_med];
Sm.ip_impact_first_iterate_closure = [first_it.closure_max];
% the mechanism at the selected state
if ~isfield(Dsel, 'vertex')
    Sm.answer1_mechanism_reproduced = false; Sm.answer1_detail = Dsel.limiting;
    Sm.answer2_centering_changed_direction = false; Sm.answer2_detail = 'no limiting iterate';
else
V = Dsel.vertex; R = Dsel.restrictions; C = Dsel.centered; I = Dsel.ipc_same_geometry;
mech = strcmp(R.binding, 'slack positivity') && Dsel.iterate.tau_hat > 0.5 && V.s_lambda_over_mu <= 0.2 ...
    && -V.predicted_motion_over_gap >= 2 && Dsel.iterate.accepted_alpha_eff <= P.short_step_threshold;
Sm.answer1_mechanism_reproduced = mech;
Sm.answer1_detail = sprintf(['frame %d, iterate k=%d (tau_hat %.3g): binding restriction %s at vertex %d (%s); s*lambda/mu = %.3g, ' ...
    'lambda/s = %.3g of the central mu/s^2; predicted normal motion %.3g x its gap; accepted fraction %.4g'], ...
    Dsel.frame, Dsel.iterate.k, Dsel.iterate.tau_hat, R.binding, V.vertex, V.source, V.s_lambda_over_mu, ...
    V.stiffness_ratio_lambda_s_to_central, V.predicted_motion_over_gap, Dsel.iterate.accepted_alpha_eff);
changed = C.limits.slack >= 2 * R.slack || C.limiting_vertex ~= V.vertex || -C.same_vertex_predicted_motion_over_gap < 1;
Sm.answer2_centering_changed_direction = changed;
Sm.answer2_detail = sprintf(['carried: slack bound %.4g (vertex %d), closure %.3g x gap; centered (lambda = mu/s): slack bound %.4g ' ...
    '(vertex %s), same vertex motion %.3g x gap, max closure %.3g x gap, direction change %.3g (rel.); IPC same geometry: ' ...
    'motion %.3g x gap, accepted %.4g'], R.slack, R.slack_vertex, -V.predicted_motion_over_gap, C.limits.slack, ...
    num2str(C.limiting_vertex), C.same_vertex_predicted_motion_over_gap, C.max_closure_over_gap, C.direction_rel_change, ...
    I.same_vertex_predicted_motion_over_gap, I.accepted_alpha);
end
ci = Sm.counts.impact; ct = Sm.counts.total;
Sm.answer3_impact_ip_more_main = ci.ip_main > ci.ipc_main;
Sm.answer3_impact_ip_more_total = ci.ip_total > ci.ipc_total;
Sm.answer3_total_ip_more_main = ct.ip_main > ct.ipc_main;
Sm.answer3_total_ip_more_total = ct.ip_total > ct.ipc_total;
Sm.failed_frames_ip = find(~[Rip.frames.qualified]);
Sm.failed_frames_ipc = find(~[Ripc.frames.qualified]);
Sm.max_rd_final_ip = max([Rip.frames.rd_final]); Sm.max_rd_final_ipc = max([Ripc.frames.rd_final]);
Sm.max_corr_final_ip = max([Rip.frames.corr_final_tol]); Sm.max_corr_final_ipc = max([Ripc.frames.corr_final_tol]);
Sm.min_gap_ip_dhat = min([Rip.frames.min_gap_dhat]); Sm.min_gap_ipc_dhat = min([Ripc.frames.min_gap_dhat]);
Sm.min_det_rel_ip = min([Rip.frames.min_det_rel]); Sm.min_det_rel_ipc = min([Ripc.frames.min_det_rel]);
Sm.max_comp_dev_ip = max([Rip.frames.comp_dev]); Sm.max_slack_consistency_ip_dhat = max([Rip.frames.slack_consistency_dhat]);
Sm.max_unresolved_ip = max([Rip.frames.unresolved]);
% event counters over the whole run
fi = Rip.frames; fc = Ripc.frames;
Sm.ip_events = struct('early_scans', sum([fi.early_scans]), 'early_rows', sum([fi.early_rows]), ...
    'early_fallbacks', sum([fi.early_fallbacks]), 'outer_clipped', sum([fi.outer_clipped]), 'split', sum([fi.split]), ...
    'dual_limited', sum([fi.dual_limited]), 'valid_limited', sum([fi.valid_limited]), 'merit_backtracks', sum([fi.backtracks]), ...
    'cert_rejects', sum([fi.cert_rejects]), 'terminal_continue', sum([fi.terminal_continue]), ...
    'redirect_would_arm', sum([fi.redirect_would_arm]), 'discarded_frames', find([fi.discarded] > 0), ...
    'rows_carried_total', sum([fi.hist_reused]), 'rows_central_at_start_total', sum([fi.hist_fresh]));
Sm.ipc_events = struct('matched_rejects', sum([fc.matched_rejects]), 'ls_limit', sum([fc.ls_limit]), ...
    'cert_rejects', sum([fc.cert_rejects]));
% impact systems by barrier level and step
imp = it(inimp);
Sm.ip_impact_systems_tau1 = sum([imp.tau_hat] > 0.5); Sm.ip_impact_systems_final = sum([imp.tau_hat] <= 0.5);
Sm.ip_impact_systems_full_step = sum([imp.alpha_eff] >= 1);
Sm.ip_impact_systems_short = sum([imp.alpha_eff] < 1 & [imp.alpha_eff] > 0);
Sm.ip_impact_systems_terminal_stop = sum([imp.alpha_eff] == 0);
d = [fi(ph.impact).main] - [fc(ph.impact).main];
Sm.impact_main_difference_per_frame = d;
Sm.impact_frames_ip_more = ph.impact(d > 0); Sm.impact_frames_ipc_more = ph.impact(d < 0);
short = ph.impact([fi(ph.impact).min_alpha_eff] <= P.short_step_threshold);
Sm.short_step_frames = short;
Sm.short_step_frames_ip_main = sum([fi(short).main]); Sm.short_step_frames_ipc_main = sum([fc(short).main]);
Sm.diagnostic_solves = Dsel.diagnostic_solves + field_or(Dmin, 'diagnostic_solves', 0);
Sm.replay_solves = Dsel.replay_solves + field_or(Dmin, 'replay_solves', 0);
Sm.wall_s = struct('ip', sum([Rip.frames.wall_s]), 'ipc', sum([Ripc.frames.wall_s]));
end

function write_summary_md(f, Sm, PT, Dsel, Dmin, Rip, Ripc, plots, SW)
fid = fopen(f, 'w');
p = @(varargin) fprintf(fid, varargin{:});
p('# Auto-generated summary (run_experiment.m)\n\n');
p('Phases: %s. First contact frame %d.\n\n', Sm.phases.rule, Sm.phases.first_contact);
p('Selected frame: %d (%s).\n\n', Sm.selection.frame, Sm.selection.rule);
p('## Newton solves by phase\n\n');
p('| phase | frames | solver | main | verification | discarded | fallback | total full-size | failed frames |\n|---|---|---|---|---|---|---|---|---|\n');
for q = 1:height(PT)
    p('| %s | %s | %s | %d | %d | %d | %d | %d | %d %s |\n', PT.phase{q}, PT.frames{q}, PT.solver{q}, PT.main_systems(q), ...
        PT.verification_systems(q), PT.discarded_systems(q), PT.fallback_systems(q), PT.total_full_size_solves(q), ...
        PT.failed_frames(q), PT.failed_frame_list{q});
end
p('\n## Impact frames, per frame\n\n| frame | IP main | IP ver | IP min accepted fraction | IPC main | IPC min accepted fraction |\n|---|---|---|---|---|---|\n');
for fr = Sm.phases.impact
    a = Rip.frames(fr); b = Ripc.frames(fr);
    p('| %d | %d | %d | %.4g | %d | %.4g |\n', fr, a.main, a.verification, a.min_alpha_eff, b.main, b.min_alpha_eff);
end
p('\n## Accuracy at the returned states\n\n');
p('- max ||r_d||_inf: IP %.4g F_w, IPC %.4g F_w (target 5)\n', Sm.max_rd_final_ip, Sm.max_rd_final_ipc);
p('- max returned correction / tol: IP %.4g (verification system), IPC %.4g (last main system)\n', Sm.max_corr_final_ip, Sm.max_corr_final_ipc);
p('- min vertex gap: IP %.4g d_hat, IPC %.4g d_hat; min det F (rel. rest): IP %.4g, IPC %.4g\n', Sm.min_gap_ip_dhat, Sm.min_gap_ipc_dhat, Sm.min_det_rel_ip, Sm.min_det_rel_ipc);
p('- IP: max |s lambda/tau - 1| %.4g, max |g - s| %.4g d_hat, max unresolved contacts %d\n', Sm.max_comp_dev_ip, Sm.max_slack_consistency_ip_dhat, Sm.max_unresolved_ip);
p('- failed frames: IP [%s], IPC [%s]\n', num2str(Sm.failed_frames_ip), num2str(Sm.failed_frames_ipc));
for D = {Dsel, Dmin}
    d = D{1}; if isempty(d), continue; end
    p('\n## Same-state replay of frame %d (%s)\n\n', d.frame, d.label);
    p('- IP replay %d + %d (trajectory %d + %d, same iterations %d, max diff alpha %.3g); IPC replay %d main (own trajectory %d); IP central start %d + %d\n', ...
        d.ip_replay.main, d.ip_replay.verification, d.ip_replay.trajectory_main, d.ip_replay.trajectory_verification, ...
        d.ip_replay.same_iteration_count, d.ip_replay.max_abs_diff_alpha, d.ipc_replay.main, d.ipc_replay.trajectory_main_on_own_trajectory, ...
        d.ip_central_start_replay.main, d.ip_central_start_replay.verification);
    if ~isfield(d, 'vertex'), p('- %s\n', d.limiting); continue; end
    V = d.vertex; R = d.restrictions; C = d.centered; I = d.ipc_same_geometry; A = d.assembly;
    p('- first limiting iterate k=%d, tau_hat %.3g, %d rows; accepted %.4g (%s)\n', d.iterate.k, d.iterate.tau_hat, d.iterate.rows, d.iterate.accepted_alpha_eff, d.iterate.limiter_recorded);
    p('- vertex %d (%s): g = %.4g d_hat, s = %.4g d_hat, lambda = %.4g F_w, mu = %.4g F_w L, s lambda/mu = %.4g, lambda/s = %.4g F_w/d_hat (central mu/s^2 = %.4g), predicted motion %.4g x gap\n', ...
        V.vertex, V.source, V.g_dhat, V.s_dhat, V.lambda_Fw, V.mu_FwL, V.s_lambda_over_mu, V.lambda_over_s_dhat_over_Fw, V.central_stiffness_mu_over_s2_dhat_over_Fw, V.predicted_motion_over_gap);
    p('- restrictions: slack %.4g (v%d), multiplier %.4g (v%s), element %.4g, ground exact crossing %.4g (v%s), merit accepted %.4g after %d backtracks; binding: %s\n', ...
        R.slack, R.slack_vertex, R.multiplier, num2str(R.multiplier_vertex), R.element, R.ground_exact, num2str(R.ground_vertex), R.merit_accepted, R.merit_backtracks, R.binding);
    p('- centered lambda = mu/s: slack %.4g (v%s), same-vertex motion %.4g x gap, max closure %.4g x gap, merit accepted %.4g, binding %s\n', ...
        C.limits.slack, num2str(C.limiting_vertex), C.same_vertex_predicted_motion_over_gap, C.max_closure_over_gap, C.limits.merit_accepted, C.limits.binding);
    p('- IPC same geometry (%d active, same pairs %d): same-vertex motion %.4g x gap, CCD %.4g (v%s), accepted %.4g\n', I.active, I.same_pairs_as_ip_rows, I.same_vertex_predicted_motion_over_gap, I.ccd, num2str(I.ccd_vertex), I.accepted_alpha);
    p('- assembly: |g - s| %.3g d_hat, stiffness rel err %.3g, rhs mu/s rel err %.3g, off-diag %.3g, KKT dx/ds/dlambda rel err %.3g / %.3g / %.3g; r_d with actual lambda %.4g F_w vs with mu/s %.4g F_w\n', ...
        A.max_abs_g_minus_s_dhat, A.stiffness_rel_err, A.rhs_effective_force_rel_err, A.offdiag_contact_terms, A.kkt_dx_rel_err, A.kkt_ds_rel_err, A.kkt_dlambda_rel_err, A.rd_inf_actual_lambda_Fw, A.rd_inf_effective_force_Fw);
    p('- diagnostic solves %d, replay solves %d (not simulation work)\n', d.diagnostic_solves, d.replay_solves);
end
p('\n## Impact sweep: every impact frame replayed from IP''s frame-start state\n\n');
p('| frame | IP main+ver | IP central start | IPC from IP state | IPC own | k | tau_hat | vertex | source | s lam/mu | motion/gap carried | slack bound carried | binding | motion/gap centered | slack bound centered | motion/gap IPC |\n');
p('|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|\n');
for q = 1:height(SW)
    p('| %d | %d+%d | %d+%d | %d | %d | %g | %.3g | %g | %s | %.3g | %.3g | %.3g | %s | %.3g | %.3g | %.3g |\n', SW.frame(q), ...
        SW.ip_main(q), SW.ip_ver(q), SW.ip_central_start_main(q), SW.ip_central_start_ver(q), SW.ipc_from_ip_state_main(q), ...
        SW.ipc_own_trajectory_main(q), SW.k(q), SW.tau_hat(q), SW.vertex(q), SW.source(q), SW.s_lambda_over_mu(q), ...
        SW.motion_over_gap_carried(q), SW.slack_bound_carried(q), SW.binding(q), SW.motion_over_gap_centered(q), ...
        SW.slack_bound_centered(q), SW.motion_over_gap_ipc(q));
end
z = Sm.impact_sweep;
p('\nTotals over the %d impact frames: IP %d + %d, IP central start %d + %d, IPC from IP''s states %d, IPC on its own trajectory %d. Diagnostic solves %d, replay solves %d.\n', ...
    z.frames, z.ip_main, z.ip_ver, z.ip_central_start_main, z.ip_central_start_ver, z.ipc_from_ip_state_main, z.ipc_own_trajectory_main, z.diagnostic_solves, z.replay_solves);
p('\n## Answers (rule-based)\n\n');
p('1. Short-step mechanism reproduced: %d. %s\n', Sm.answer1_mechanism_reproduced, Sm.answer1_detail);
p('2. Centering changed the limiting direction: %d. %s\n', Sm.answer2_centering_changed_direction, Sm.answer2_detail);
ci = Sm.counts.impact; ct = Sm.counts.total;
p('3. Impact: IP %d main + %d ver = %d vs IPC %d. Full run: IP %d main + %d ver = %d vs IPC %d.\n', ci.ip_main, ci.ip_verification, ci.ip_total, ci.ipc_main, ct.ip_main, ct.ip_verification, ct.ip_total, ct.ipc_main);
p('\nPlots: %s\n', strjoin(plots, ', '));
fclose(fid);
end

function v = field_or(s, f, d)
if isempty(s) || ~isfield(s, f), v = d; else, v = s.(f); end
end
