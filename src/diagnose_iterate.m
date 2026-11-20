function L = diagnose_iterate(S, P, it)
%DIAGNOSE_ITERATE  Step-limit analysis of one captured hard-IP iterate (it.cap from ip_solve_frame(..., capture=true)).
% The step is re-formed with the solver's own shared rules (ip_step_proposal, ip_line_search with ip_armijo and its
% roundoff rule, ip_step_limiter, the all-vertex ground scan after the accept), so the reported limiter is the one that
% actually set the step, whatever it is. Contact fields (g, s, lambda, s lambda/mu, lambda/s, predicted motion) are
% filled only when a contact row (or, for a ground fallback, a vertex) sets the step; otherwise they are empty.
% Also: the assembly identities, an unreduced KKT solve, one diagnostic direction with lambda_i := mu/s_i at the same
% positions, rows, slacks and mu, and IPC's own direction and step at the same geometry.
% Diagnostic solves are counted in L.diagnostic_solves (they are not simulation work).
c = it.cap;
assert(~isempty(c), 'diagnose_iterate: the iterate has no capture (run ip_solve_frame with capture = true)');
Fw = S.F_w; dhat = S.d_hat; tol = S.tol_disp; gy = S.ground_y;
x = c.x; xt = c.xt; rv = c.rv; s = c.s; lam = c.lam; tau = c.tau; g = c.g; mu = tau;
y = x(S.ydof) - gy;
L = struct('found', true, 'k', it.k, 'solve', it.solve, 'tau_hat', it.tau_hat, 'rows', numel(rv), 'levels', c.levels, ...
    'diagnostic_solves', 0);
L.recorded = struct('alpha', it.alpha, 'alpha_eff', it.alpha_eff, 'fallback_t', it.fallback_t, 'limiter', it.limiter, ...
    'lim_vertex', it.lim_vertex, 'alpha_slack', it.alpha_slack, 'alpha_mult', it.alpha_mult, 'alpha_elem', it.alpha_elem, ...
    'alpha_proposal', it.alpha_proposal, 'trials', it.trials);

% ---- the direction, recomputed
[E, gE, H0] = mech_eval(S, x, xt, 3);
JT = @(r) accum_y(S, rv, r);
w = lam ./ s;
K = H0 + sparse(2*rv, 2*rv, w, S.ndof, S.ndof);
b = gE - JT(mu ./ s) + JT(w .* (g - s));
dx = shared_solve(K, -b); L.diagnostic_solves = L.diagnostic_solves + 1;
L.recomputed_direction_max_abs_diff = max(abs(dx - c.dx));
dx = c.dx;   % analyze the solver's own direction
ds = dx(2*rv) + (g - s);
dlam = mu ./ s - lam - w .* ds;

% ---- assembly identities and the unreduced Newton system
A = struct();
A.max_abs_g_minus_s_dhat = maxabs(g - s) / dhat;
A.s_equals_g = isempty(rv) || all(g == s);
dK = K - H0;
A.stiffness_rel_err = relmax(full(diag(dK(2*rv, 2*rv))), lam ./ s);
dK(sub2ind(size(dK), 2*rv, 2*rv)) = 0;
A.offdiag_contact_terms = full(max(abs(dK(:))));
db = b - gE;
A.rhs_general_rel_err = relmax(db(2*rv), -mu ./ s + w .* (g - s));       % holds for any s
if A.s_equals_g
    A.rhs_effective_force_rel_err = relmax(-db(2*rv), mu ./ s);          % the s = g form: contact force mu/s
else
    A.rhs_effective_force_rel_err = [];                                   % not applicable: s ~= g at this iterate
end
rd = gE - JT(lam);
A.rd_inf_actual_lambda_Fw = max(abs(rd)) / Fw;
A.rd_inf_effective_force_Fw = max(abs(gE - JT(mu ./ s))) / Fw;
m = numel(rv); n = S.ndof;
Jm = sparse(1:m, 2*rv, 1, m, n);
KKT = [H0, sparse(n, m), -Jm'; Jm, -speye(m), sparse(m, m); sparse(m, n), spdiags(lam, 0, m, m), spdiags(s, 0, m, m)];
z = KKT \ (-[rd; g - s; s .* lam - mu]); L.diagnostic_solves = L.diagnostic_solves + 1;
A.kkt_dx_rel_err = norm(z(1:n) - dx) / norm(dx);
A.kkt_ds_rel_err = norm(z(n+1:n+m) - ds) / max(norm(ds), eps);
A.kkt_dlambda_rel_err = norm(z(n+m+1:end) - dlam) / max(norm(dlam), eps);
L.assembly = A;

% ---- the carried-multiplier step, re-formed with the solver's rules
Rc = analyze_direction(S, P, x, xt, rv, s, lam, tau, g, y, dx, ds, dlam, c.nu, c.slope, E, c.levels);
Rc.agrees_with_solver = strcmp(Rc.limiter, it.limiter) && same_num(Rc.accepted_alpha_eff, it.alpha_eff);
L.carried = Rc;
L.limiter = Rc.limiter;
L.vertex = vertex_fields(S, P, Rc, rv, s, lam, mu, g, y, dx, ds, dlam, gE, H0, c.src);

% ---- diagnostic direction: lambda_i := mu / s_i, everything else identical
lam_c = mu ./ s; w_c = lam_c ./ s;
K_c = H0 + sparse(2*rv, 2*rv, w_c, S.ndof, S.ndof);
b_c = gE - JT(mu ./ s) + JT(w_c .* (g - s));
dx_c = shared_solve(K_c, -b_c); L.diagnostic_solves = L.diagnostic_solves + 1;
ds_c = dx_c(2*rv) + (g - s);
dlam_c = mu ./ s - lam_c - w_c .* ds_c;
base_c = gE' * dx_c - tau * sum(ds_c ./ s);
[nu_c, slope_c] = ip_merit_slope(base_c, sum(abs(g - s)), maxabs(lam_c + dlam_c), c.nu_prev, P.ip);
Rz = analyze_direction(S, P, x, xt, rv, s, lam_c, tau, g, y, dx_c, ds_c, dlam_c, nu_c, slope_c, E, c.levels);
C = struct('limits', Rz, 'limiter', Rz.limiter, 'limiting_vertex', Rz.limiting_vertex, ...
    'direction_rel_change', norm(dx_c - dx) / norm(dx), 'max_closure_over_gap', maxval(-dx_c(2*rv) ./ g), ...
    'rd_inf_with_centered_lambda_Fw', max(abs(gE - JT(lam_c))) / Fw, 'same_vertex_predicted_motion_over_gap', [], ...
    'same_vertex_Jdx_dhat', [], 'same_vertex_stiffness_lambda_over_s_dhat_over_Fw', []);
if L.vertex.applicable && L.vertex.has_row
    j = L.vertex.row;
    C.same_vertex_predicted_motion_over_gap = dx_c(2*rv(j)) / g(j);
    C.same_vertex_Jdx_dhat = dx_c(2*rv(j)) / dhat;
    C.same_vertex_stiffness_lambda_over_s_dhat_over_Fw = w_c(j) * dhat / Fw;
end
L.centered = C;

% ---- IPC at the same geometry (its own active set, force and stiffness at x)
[Em, gm, Hm] = mech_eval(S, x, xt, 3);
[Eb, gb, Hb, Cb] = ipc_contact_eval(S, x, 3);
dx_i = shared_solve(Hm + Hb, -(gm + gb)); L.diagnostic_solves = L.diagnostic_solves + 1;
I = struct();
I.active = numel(Cb.act); I.same_pairs_as_ip_rows = isequal(sort(Cb.act), sort(rv));
I.rd_inf_Fw = max(abs(gm + gb)) / Fw; I.dx_tol = max(abs(dx_i)) / tol;
dyi = dx_i(S.ydof); idx = find(dyi < 0);
[tm, qq] = min(y(idx) ./ (-dyi(idx)) * (1 - P.ipc.ccd_eta));
if isempty(tm) || tm >= 1.1, I.ccd = 1; I.ccd_vertex = NaN; else, I.ccd = min(1, tm); I.ccd_vertex = idx(qq); end
fe = elem_check(S, x, dx_i, I.ccd);
I.elem = Inf; if fe <= I.ccd, I.elem = P.valid_eta * fe; end
a = min(I.ccd, I.elem); E0 = Em + Eb; I.trials = zeros(0, 3);
for t = 1:P.ls_max
    xtr = x + a * dx_i;
    [~, np] = elem_check(S, xtr);
    if np > 0, I.trials = [I.trials; a, NaN, 0]; a = a / 2; continue; end
    Et = mech_eval(S, xtr, xt, 1) + ipc_contact_eval(S, xtr, 1);
    I.trials = [I.trials; a, Et - E0, Et <= E0];
    if Et <= E0, break; end
    a = a / 2;
end
I.accepted_alpha = a;
I.max_closure_over_gap_ip_rows = maxval(-dx_i(2*rv) ./ g);
I.same_vertex_predicted_motion_over_gap = []; I.same_vertex_Jdx_dhat = [];
I.same_vertex_stiffness_dhat_over_Fw = []; I.same_vertex_force_Fw = [];
if L.vertex.applicable
    lv = L.vertex.vertex;
    I.same_vertex_predicted_motion_over_gap = dx_i(2*lv) / y(lv);
    I.same_vertex_Jdx_dhat = dx_i(2*lv) / dhat;
    ia = find(Cb.act == lv);
    if ~isempty(ia), I.same_vertex_stiffness_dhat_over_Fw = Cb.kn(ia) * dhat / Fw; I.same_vertex_force_Fw = Cb.fn(ia) / Fw; end
end
L.ipc_same_geometry = I;

% ---- every row: carried, centered, IPC
L.rows_table = table(rv, g / dhat, s / dhat, lam / Fw, s .* lam / mu, w * dhat / Fw, (mu ./ s.^2) * dhat / Fw, ...
    dx(2*rv) ./ g, dx_c(2*rv) ./ g, dx_i(2*rv) ./ g, ftb_each(s, ds, P.ip.ftb_eta), ftb_each(s, ds_c, P.ip.ftb_eta), ...
    ftb_each(lam, dlam, P.ip.ftb_eta), arrayfun(@(v) ipc_kn(S, y(v)), rv) * dhat / Fw, ...
    reshape(arrayfun(@(q2) string(src_name(c.src(q2))), 1:numel(rv)), [], 1), ...
    'VariableNames', {'vertex', 'g_dhat', 's_dhat', 'lambda_Fw', 's_lambda_over_mu', 'lambda_over_s_dhat_Fw', ...
    'mu_over_s2_dhat_Fw', 'motion_over_gap_carried', 'motion_over_gap_centered', 'motion_over_gap_ipc', ...
    'slack_bound_carried', 'slack_bound_centered', 'multiplier_bound_carried', 'ipc_kn_dhat_Fw', 'source'});
end

% ======================================================================================== the step, as the solver forms it
function R = analyze_direction(S, P, x, xt, rv, s, lam, tau, g, y, dx, ds, dlam, nu, slope, E, levels)
SP = ip_step_proposal(S, P, x, s, lam, ds, dlam, dx, levels);
phi0 = E - tau * sum(log(s / S.d_hat)) + nu * sum(abs(g - s));
LS = ip_line_search(S, P, x, xt, rv, s, lam, dx, ds, dlam, tau, nu, slope, phi0, SP.alpha_proposal, SP.split_now, SP.a_l);
[name, row] = ip_step_limiter(SP, LS);
R = struct();
R.slack = SP.a_s; R.decreasing_slacks = sum(ds < 0);
R.multiplier = SP.a_l; R.decreasing_multipliers = sum(dlam < 0);
% a bound names a vertex only when it restricts the step (fraction < 1)
R.slack_vertex = NaN; if SP.a_s < 1, R.slack_vertex = vid(rv, SP.js); end
R.multiplier_vertex = NaN; if SP.a_l < 1, R.multiplier_vertex = vid(rv, SP.jl); end
R.split = SP.split_now; R.update_B_multiplier_limited = SP.dual_lim && ~SP.split_now;
R.element = SP.a_elem; R.element_triangle = SP.elem_tri; R.proposal = SP.alpha_proposal;
dy = dx(S.ydof); idx = find(dy < 0);
[tg, qg] = min(y(idx) ./ (-dy(idx)));
if isempty(tg), R.ground_exact = Inf; R.ground_vertex = NaN; else, R.ground_exact = tg; R.ground_vertex = idx(qg); end
R.ground_vertex_has_row = ismember(R.ground_vertex, rv);
R.line_search = struct('accepted', LS.accepted, 'alpha', LS.alpha, 'trials', LS.trials, 'element_rejects', LS.cert, ...
    'by_precision', LS.by_precision, 'log', LS.log, ...
    'log_columns', 'alpha, phi - phi0, armijo_c alpha slope, accepted, accepted by the roundoff rule, element-rejected');
R.merit_accepted = NaN; if LS.accepted, R.merit_accepted = LS.alpha; end
R.merit_backtracks = LS.trials - 1;
R.limiter = name; R.limiting_row = row; R.limiting_vertex = vid(rv, row);
R.fallback_t = NaN; R.scan_all_positive = NaN; R.min_gap_after_dhat = NaN;
R.accepted_alpha = 0; R.accepted_alpha_eff = 0;
if LS.accepted
    R.accepted_alpha = LS.alpha; R.accepted_alpha_eff = LS.alpha;
    yacc = LS.x(S.ydof) - S.ground_y;
    R.scan_all_positive = all(yacc > 0); R.min_gap_after_dhat = min(yacc) / S.d_hat;
    if ~R.scan_all_positive
        % the early-scan fallback: halving from the step's start until every vertex gap is positive
        tgf = 1; h = 0; seg = LS.alpha * dy;
        while ~all(y + tgf * seg > 0) && h < P.ip.ground_halvings, tgf = tgf / 2; h = h + 1; end
        R.fallback_t = tgf; R.accepted_alpha_eff = LS.alpha * tgf;
        [~, R.limiting_vertex] = min(yacc); R.limiting_row = find(rv == R.limiting_vertex, 1);
        if isempty(R.limiting_row), R.limiting_row = 0; end
        R.limiter = [name ' + ground fallback'];
    end
end
R.binding = R.limiter;
end

function V = vertex_fields(S, P, R, rv, s, lam, mu, g, y, dx, ds, dlam, gE, H0, src)
% contact quantities of the vertex whose restriction set the step; empty when no contact row or vertex did
Fw = S.F_w; dhat = S.d_hat;
V = struct('applicable', false, 'reason', '', 'vertex', [], 'has_row', false, 'row', [], 'source', '', ...
    'is_bottom_vertex', [], 'g_dhat', [], 's_dhat', [], 'lambda_Fw', [], 'mu_FwL', mu / (Fw * S.L), ...
    's_lambda_over_mu', [], 'lambda_over_s_dhat_over_Fw', [], 'lambda_over_s_over_mv', [], ...
    'central_stiffness_mu_over_s2_dhat_over_Fw', [], 'stiffness_ratio_lambda_s_to_central', [], ...
    'effective_force_mu_over_s_Fw', [], 'ipc_force_same_gap_Fw', [], 'ipc_stiffness_same_gap_dhat_over_Fw', [], ...
    'Jdx_dhat', [], 'predicted_motion_over_gap', [], 'ds_dhat', [], 'dlambda_Fw', [], 'raw_slack_crossing', [], ...
    'slack_bound', [], 'multiplier_bound', [], 'noncontact_force_Fw', [], 'noncontact_stiffness_over_mv', []);
if isnan(R.limiting_vertex)
    V.reason = sprintf('no contact row or vertex sets this step (limiter: %s)', R.limiter);
    return;
end
lv = R.limiting_vertex;
V.applicable = true; V.vertex = lv; V.is_bottom_vertex = ismember(lv, S.bottom);
V.g_dhat = y(lv) / dhat;
V.Jdx_dhat = dx(2*lv) / dhat; V.predicted_motion_over_gap = dx(2*lv) / y(lv);
V.ipc_force_same_gap_Fw = ipc_fn(S, y(lv)) / Fw; V.ipc_stiffness_same_gap_dhat_over_Fw = ipc_kn(S, y(lv)) * dhat / Fw;
V.noncontact_force_Fw = -gE(2*lv) / Fw; V.noncontact_stiffness_over_mv = full(H0(2*lv, 2*lv)) / S.m(lv);
j = R.limiting_row;
if j == 0
    V.reason = 'the limiting vertex has no contact row (ground fallback on an untracked vertex)';
    return;
end
V.has_row = true; V.row = j;
V.s_dhat = s(j) / dhat; V.lambda_Fw = lam(j) / Fw;
V.s_lambda_over_mu = s(j) * lam(j) / mu;
V.lambda_over_s_dhat_over_Fw = lam(j) / s(j) * dhat / Fw;
V.lambda_over_s_over_mv = lam(j) / s(j) / S.m(lv);
V.central_stiffness_mu_over_s2_dhat_over_Fw = mu / s(j)^2 * dhat / Fw;
V.stiffness_ratio_lambda_s_to_central = (lam(j) / s(j)) / (mu / s(j)^2);
V.effective_force_mu_over_s_Fw = mu / s(j) / Fw;
V.ds_dhat = ds(j) / dhat; V.dlambda_Fw = dlam(j) / Fw;
if ds(j) < 0
    V.raw_slack_crossing = s(j) / (-ds(j)); V.slack_bound = P.ip.ftb_eta * s(j) / (-ds(j));
end
if dlam(j) < 0
    V.multiplier_bound = P.ip.ftb_eta * lam(j) / (-dlam(j));
end
V.source = src_name(src(j));
end

% ======================================================================================== helpers
function v = accum_y(S, rv, r)
v = zeros(S.ndof, 1); v(2*rv) = r;
end

function e = relmax(a, b)
if isempty(a), e = []; return; end
e = max(abs(a - b) ./ max(abs(b), realmin));
end

function m = maxabs(v)
if isempty(v), m = 0; else, m = max(abs(v)); end
end

function m = maxval(v)
if isempty(v), m = []; else, m = max(v); end
end

function v = vid(rv, j)
if j > 0, v = rv(j); else, v = NaN; end
end

function a = ftb_each(v, dv, eta)
a = inf(size(v)); q = dv < 0;
a(q) = eta * v(q) ./ (-dv(q));
end

function tf = same_num(a, b)
tf = (isnan(a) && isnan(b)) || a == b;
end

function f = ipc_fn(S, gap)
if gap >= S.d_hat || gap <= 0, f = 0; return; end
[~, dB] = ipc_barrier(gap^2, S.d_hat^2);
f = -2 * S.kappa_dt2 * dB * gap;
end

function k = ipc_kn(S, gap)
if gap >= S.d_hat || gap <= 0, k = 0; return; end
[~, dB, ddB] = ipc_barrier(gap^2, S.d_hat^2);
k = S.kappa_dt2 * (4 * gap^2 * ddB + 2 * dB);
end

function n = src_name(code)
names = {'central (frame start)', 'carried (history)', 'early insertion (central)', 'outer refresh (central)'};
n = names{code + 1};
end
