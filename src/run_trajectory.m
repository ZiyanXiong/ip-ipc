function R = run_trajectory(S, P, solver, frames, opts)
%RUN_TRAJECTORY  One solver on its own trajectory from the scene start.
% Implicit Euler: xt = x_n + dt v_n + dt^2 g (shared inertial prediction), v_{n+1} = (x_{n+1} - x_n) / dt.
%
% Failure handling. A frame qualifies only if its solver returns status 'converged' (IPC: the matched stop; IP: the
% terminal update and the fresh verification). The first frame that does not qualify ENDS the trajectory:
%   - its diagnostics are kept in R.failure (frame record, every iteration, the unqualified returned state);
%   - the state is not advanced: R.x_end, R.v_end and (IP) R.hist_end are those of the last qualified frame, and
%     physical time stops there (R.time_end = R.last_qualified_frame * dt);
%   - R.complete = false and R.first_failure = that frame. Such a trajectory is excluded from performance
%     comparisons (compare_performance.m).
% opts.fault (tests only): struct('frame', f, 'cap_main', c) replaces the main-system budget of frame f, to exercise
% the failure path through the solver's own budget stop. Normal runs never set it.
if nargin < 4 || isempty(frames), frames = P.frames; end
if nargin < 5 || isempty(opts), opts = struct(); end
fault = [];
if isfield(opts, 'fault'), fault = opts.fault; end
x = S.x0; v = zeros(S.ndof, 1); hist = [];
fr = []; its = []; starts = {};
failure = []; last_q = 0; first_failure = NaN;
for f = 1:frames
    xt = x + S.dt * v + S.dt2 * S.gvec;
    starts{end+1, 1} = struct('x', x, 'v', v, 'xt', xt, 'hist', hist); %#ok<AGROW>
    Pf = P;
    injected = ~isempty(fault) && fault.frame == f;
    if injected, Pf.cap_main = fault.cap_main; end
    t0 = tic;
    switch solver
        case 'ip'
            [xn, F, it, hist_out] = ip_solve_frame(S, Pf, x, xt, f, hist, false);
        case 'ipc'
            [xn, F, it] = ipc_solve_frame(S, Pf, x, xt, f, false);
            hist_out = [];
        otherwise
            error('run_trajectory:solver', 'unknown solver %s', solver);
    end
    F.wall = toc(t0);
    F.contact = any([it.rows] > 0);
    F.qualified = strcmp(F.status, 'converged');
    F.fault_injected = injected;
    % metrics of the state this frame returned (for a failed frame: the unqualified attempt, not adopted)
    y = xn(S.ydof) - S.ground_y;
    [~, np, dets] = elem_check(S, xn);
    feasible = all(y > 0) && np == 0;
    F.min_gap_dhat = min(y) / S.d_hat;
    F.min_det_rel = min(dets ./ (2 * S.A(:)));
    if F.qualified
        v = (xn - x) / S.dt; x = xn;
        if strcmp(solver, 'ip'), hist = hist_out; end
        last_q = f;
    end
    F.state_advanced = F.qualified;
    X = reshape(x, 2, []);
    F.com_x = sum(S.m' .* X(1,:)) / sum(S.m); F.com_y = sum(S.m' .* X(2,:)) / sum(S.m);
    F.vy_com = sum(S.m .* v(S.ydof)) / sum(S.m);
    F.total = F.main + F.verification;
    rec = frame_summary(F, it);
    fr = [fr, rec]; %#ok<AGROW>
    for q = 1:numel(it), it(q).cap = []; end
    its = [its, it]; %#ok<AGROW>
    if ~F.qualified
        first_failure = f;
        failure = struct('frame', f, 'status', F.status, 'reason', F.reason, 'fault_injected', injected, ...
            'frame_record', rec, 'iters', it, 'returned_state', xn, 'returned_state_feasible', feasible, ...
            'frame_start', starts{end});
        break;
    end
end
R = struct('solver', solver, 'frames', fr, 'iters', its, 'starts', {starts}, 'x_end', x, 'v_end', v, ...
    'hist_end', {hist}, 'complete', last_q == frames, 'frames_requested', frames, 'last_qualified_frame', last_q, ...
    'first_failure', first_failure, 'time_end', last_q * S.dt, 'failure', {failure});
end

function o = frame_summary(F, it)
o = struct();
o.solver = F.solver; o.frame = F.frame; o.status = F.status; o.reason = F.reason; o.qualified = F.qualified;
o.state_advanced = F.state_advanced; o.fault_injected = F.fault_injected;
o.main = F.main; o.verification = F.verification; o.discarded = F.discarded; o.fallback = F.fallback;
o.total = F.main + F.verification; o.wall_s = F.wall; o.contact = F.contact; o.trials = F.trials;
if isempty(it), o.rows_or_active_end = NaN; else, o.rows_or_active_end = it(end).rows; end
o.min_gap_dhat = F.min_gap_dhat; o.min_det_rel = F.min_det_rel;
o.com_x = F.com_x; o.com_y = F.com_y; o.vy_com = F.vy_com;
if isempty(it), mains = it; else, mains = it(strcmp({it.kind}, 'main')); end
acc = [mains.alpha_eff]; acc = acc(isfinite(acc) & acc > 0);
if isempty(acc), o.min_alpha_eff = NaN; else, o.min_alpha_eff = min(acc); end
switch F.solver
    case 'ip'
        o.ls_limit = NaN; o.matched_rejects = NaN;
        o.levels = F.levels; o.split = F.split; o.dual_limited = F.dual_limited; o.valid_limited = F.valid_limited;
        o.early_scans = F.early_scans; o.early_rows = F.early_rows; o.early_fallbacks = F.early_fallbacks;
        o.outer_clipped = F.outer_clipped; o.backtracks = F.backtracks; o.cert_rejects = F.cert_rejects;
        o.terminal_applied = F.terminal_applied; o.terminal_continue = F.terminal_continue;
        o.hist_reused = F.hist_reused; o.hist_fresh = F.hist_fresh; o.rows_start = F.rows_start;
        o.redirect_would_arm = F.redirect_would_arm;
        v = F.ver;
        if isempty(v)
            o.rd_final = NaN; o.corr_final_tol = NaN; o.slack_consistency_dhat = NaN; o.comp_dev = NaN;
            o.unresolved = NaN; o.ver_ok = false;
        else
            o.rd_final = v.rd_inf; o.corr_final_tol = v.dx_tol; o.slack_consistency_dhat = v.slack_consistency_dhat;
            o.comp_dev = v.comp_dev; o.unresolved = v.unresolved; o.ver_ok = v.hit;
        end
    case 'ipc'
        o.ls_limit = F.ls_limit; o.matched_rejects = F.matched_rejects;
        o.levels = NaN; o.split = NaN; o.dual_limited = NaN; o.valid_limited = NaN;
        o.early_scans = NaN; o.early_rows = NaN; o.early_fallbacks = NaN; o.outer_clipped = NaN;
        o.backtracks = NaN; o.cert_rejects = F.cert_rejects;
        o.terminal_applied = NaN; o.terminal_continue = NaN; o.hist_reused = NaN; o.hist_fresh = NaN;
        if isempty(it), o.rows_start = NaN; else, o.rows_start = it(1).rows; end
        o.redirect_would_arm = NaN;
        if isempty(it)
            o.rd_final = NaN; o.corr_final_tol = NaN;
        else
            o.rd_final = it(end).rd_inf; o.corr_final_tol = it(end).dx_tol;
        end
        o.slack_consistency_dhat = NaN; o.comp_dev = NaN; o.unresolved = NaN;
        o.ver_ok = strcmp(F.status, 'converged');
end
end
