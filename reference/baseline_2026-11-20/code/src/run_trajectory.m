function R = run_trajectory(S, P, solver, frames)
%RUN_TRAJECTORY  One solver on its own trajectory from the scene start.
% Implicit Euler: xt = x_n + dt v_n + dt^2 g (shared inertial prediction), v_{n+1} = (x_{n+1} - x_n) / dt.
% A frame that does not qualify (status ~= converged) is recorded; the trajectory continues from its last iterate if
% that iterate is feasible (all gaps > 0, valid elements); otherwise from the frame start with no velocity update.
% A failed IP frame leaves no multiplier history (the next frame starts centrally).
if nargin < 4, frames = P.frames; end
x = S.x0; v = zeros(S.ndof, 1); hist = [];
fr = []; its = []; starts = cell(frames, 1);
for f = 1:frames
    xt = x + S.dt * v + S.dt2 * S.gvec;
    starts{f} = struct('x', x, 'v', v, 'xt', xt, 'hist', hist);
    t0 = tic;
    switch solver
        case 'ip'
            [xn, F, it, hist_out] = ip_solve_frame(S, P, x, xt, f, hist, false);
        case 'ipc'
            [xn, F, it] = ipc_solve_frame(S, P, x, xt, f, false);
            hist_out = [];
    end
    F.wall = toc(t0);
    F.contact = F_contact(it);
    y = xn(S.ydof) - S.ground_y;
    [~, np, dets] = elem_check(S, xn);
    feasible = all(y > 0) && np == 0;
    F.qualified = strcmp(F.status, 'converged');
    F.continued_from = 'returned state';
    if F.qualified || feasible
        v = (xn - x) / S.dt; x = xn;
        if ~F.qualified, F.continued_from = 'last iterate (unqualified)'; end
    else
        F.continued_from = 'frame start (no velocity update)';
    end
    if strcmp(solver, 'ip')
        if F.qualified, hist = hist_out; else, hist = []; end
    end
    F.min_gap_dhat = min(y) / S.d_hat;
    F.min_det_rel = min(dets ./ (2 * S.A(:)));
    X = reshape(x, 2, []);
    F.com_x = sum(S.m' .* X(1,:)) / sum(S.m); F.com_y = sum(S.m' .* X(2,:)) / sum(S.m);
    F.vy_com = sum(S.m .* v(S.ydof)) / sum(S.m);
    F.total = F.main + F.verification;
    fr = [fr, frame_summary(S, F, it)]; %#ok<AGROW>
    for q = 1:numel(it), it(q).cap = []; end
    its = [its, it]; %#ok<AGROW>
end
R = struct('solver', solver, 'frames', fr, 'iters', its, 'starts', {starts}, 'x_end', x, 'v_end', v);
end

function c = F_contact(it)
c = any([it.rows] > 0);
end

function o = frame_summary(S, F, it)
o = struct();
o.solver = F.solver; o.frame = F.frame; o.status = F.status; o.reason = F.reason; o.qualified = F.qualified;
o.continued_from = F.continued_from;
o.main = F.main; o.verification = F.verification; o.discarded = F.discarded; o.fallback = F.fallback;
o.total = F.main + F.verification; o.wall_s = F.wall; o.contact = F.contact; o.trials = F.trials;
o.rows_or_active_end = it(end).rows;
o.min_gap_dhat = F.min_gap_dhat; o.min_det_rel = F.min_det_rel;
o.com_x = F.com_x; o.com_y = F.com_y; o.vy_com = F.vy_com;
mains = it(strcmp({it.kind}, 'main'));
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
        o.rows_start = it(1).rows; o.redirect_would_arm = NaN;
        last = it(end);
        o.rd_final = last.rd_inf; o.corr_final_tol = last.dx_tol;
        o.slack_consistency_dhat = NaN; o.comp_dev = NaN; o.unresolved = NaN;
        o.ver_ok = strcmp(F.status, 'converged');
end
end
