function [x, F, its, hist_out, st] = ip_solve_frame(S, P, xn, xt, frame, hist, capture)
%IP_SOLVE_FRAME  One frame of ordinary hard IP, the round-43 'base' configuration (2D port of engine/ip_solver.cu):
% arm H (multiplier history, terminal update) + signed-gap schedule + zero-CCD ground advance + early discovery.
%
% Rows i are vertex-ground pairs with gap g_i(x) = y_i, slack s_i > 0, multiplier lambda_i > 0, barrier tau (mu):
%   g_i - s_i = 0,   s_i lambda_i = tau,   grad E_mech - J' lambda = 0.
% Newton, with ds and dlambda eliminated (the plane has zero curvature, so the member block is (lambda/s) n n'):
%   (H0 + J' diag(lambda/s) J) dx = -(grad E - J' (tau/s) + J' diag(lambda/s) (g - s))
%   ds = J dx + (g - s),   dlambda = tau/s - lambda - (lambda/s) ds
% Step: slack and multiplier fractions to the boundary (eta 0.995); update B (one fraction) except at the first
% barrier level, where the split rule moves positions and slacks by the slack fraction and the multipliers by
% min(alpha, their own fraction); element-validity bound; Armijo on phi = F - tau sum log(s~/d_hat) + nu sum|g - s~|
% with s~ = max(s + alpha ds, g(x~)), nu = 1.1 max|lambda + dlambda| (descent safeguard); halving, ls_max trials.
% No CCD. After every accepted step the gap of EVERY vertex is checked (early discovery): all positive -> untracked
% vertices inside the band (g < d_hat) become rows (s = g, lambda = tau/g, central), the point becomes the last checked
% position x_k; any nonpositive -> the ground advance from x_k (full, then halving until every gap > 0; element
% bound), every row's slack reset to its gap there (multipliers kept), refresh.
% Barrier: tau_hat = 1 at the frame start; the first decrease goes directly to the final level tol_disp / L.
% Stop: final level, |dx|_inf <= tol, |g - s|_inf <= tol, |s lambda - tau|_inf <= kappa_eps tau_final, the full
% correction keeps every vertex above ground and every element valid -> terminal update: the correction goes through
% the ordinary safeguarded update, then the fresh verification of the returned state (one verification system).
if nargin < 6, hist = []; end
if nargin < 7, capture = false; end
tol = S.tol_disp; Fw = S.F_w; L = S.L; dhat = S.d_hat; gy = S.ground_y; ydof = S.ydof; Q = P.ip;
x = xn;
tau = S.tau0; levels = 0;
F = struct('solver', 'ip', 'frame', frame, 'status', 'running', 'reason', '', 'main', 0, 'verification', 0, ...
    'discarded', 0, 'fallback', 0, 'accepted', 0, 'levels', 0, 'split', 0, 'dual_limited', 0, 'valid_limited', 0, ...
    'cert_rejects', 0, 'backtracks', 0, 'trials', 0, 'early_scans', 0, 'early_batches', 0, 'early_rows', 0, ...
    'early_fallbacks', 0, 'outer_events', 0, 'outer_clipped', 0, 'refresh_rows', 0, 'nu_safeguard', 0, ...
    'terminal_candidates', 0, 'terminal_applied', 0, 'terminal_rejected', 0, 'terminal_continue', 0, ...
    'hist_used', false, 'hist_why', '', 'hist_reused', 0, 'hist_fresh', 0, 'rows_start', 0, 'rows_end', 0, ...
    'redirect_would_arm', 0, 'ver', []);
its = []; hist_out = []; st = struct();

% ---- frame-start discovery, central start, multiplier history
y = x(ydof) - gy;
if any(y <= 0)
    F.status = 'failed'; F.reason = 'a vertex is not above the ground at the frame start'; return;
end
rv = find(y < dhat);
s = max(y(rv), S.s_floor);
lam = tau ./ s;
src = zeros(size(rv));             % 0 central, 1 reused, 2 early insertion, 3 outer refresh
born = zeros(size(rv)) - 1;
if Q.history
    if isempty(hist)
        F.hist_why = 'no previous accepted frame';
    elseif hist.frame + 1 ~= frame
        F.hist_why = 'the history is not the immediately preceding accepted frame';
    else
        F.hist_used = true;
        [tf, loc] = ismember(rv, hist.rv);
        for q = find(tf)'
            hv = hist.lam(loc(q));
            if isfinite(hv) && hv > 0
                lam(q) = hv; src(q) = 1;
            end
        end
    end
end
F.hist_reused = sum(src == 1); F.hist_fresh = numel(rv) - F.hist_reused; F.rows_start = numel(rv);
xk = x;                  % the last checked (collision-safe) outer position
so_early_skip = true;    % the start position was just discovered
so_refresh = false; x_moved = false;
last_dx = -1; nu = 0; stalled = 0; consec_negl = 0; lim_hist = {};
ver_last = [];

for k = 0:P.cap_main-1
    ev = {}; inserted = 0;
    % ---------------- iteration top: early ground discovery at the accepted point
    if x_moved
        if so_refresh
            so_refresh = false;
            [ok, nadd] = scan_or_fallback('refresh after a terminal continue');
            if ~ok, break; end
            inserted = nadd;
        elseif Q.early_discovery && ~so_early_skip
            [ok, nadd] = scan_or_fallback('early');
            if ~ok, break; end
            inserted = nadd;
        end
    end
    so_early_skip = false; x_moved = false;

    % ---------------- residuals at this point
    y = x(ydof) - gy; g = y(rv);
    [E, gE, H0] = mech_eval(S, x, xt, 3);
    rd = gE - JT(lam);
    rd_inf = max(abs(rd));
    rp = g - s;
    rp_inf = maxabs(rp);
    cent = maxabs(s .* lam - tau);
    % ---------------- residual-driven barrier continuation (first decrease direct to the final level)
    while should_decrease(tau, rd_inf, rp_inf, cent, last_dx)
        if Q.tau_direct && levels == 0
            tn = S.tau_final;
        else
            th = tau / (Fw * L);
            tn = max(S.tau_final, min(Q.kappa_tau * th, th^Q.theta_tau) * Fw * L);
        end
        ev{end+1} = sprintf('tau_hat %.4g -> %.4g', tau / (Fw*L), tn / (Fw*L)); %#ok<AGROW>
        levels = levels + 1; tau = tn; last_dx = -1;
        cent = maxabs(s .* lam - tau);
    end
    % ---------------- condensed Newton system (one main system)
    w = lam ./ s;
    K = H0 + sparse(2*rv, 2*rv, w, S.ndof, S.ndof);
    b = gE - JT(tau ./ s) + JT(w .* rp);
    dx = shared_solve(K, -b);
    F.main = F.main + 1;
    dy = dx(ydof);
    Jdx = dy(rv);
    ds = Jdx + rp;
    dlam = tau ./ s - lam - w .* ds;
    dx_inf = max(abs(dx));
    % merit slope and weight
    gF_dx = gE' * dx;
    base = gF_dx - tau * sum(ds ./ s);
    rp1 = sum(abs(rp));
    lam_plus = maxabs(lam + dlam);
    est = Q.nu_rho * lam_plus;
    if isfinite(est), nu = est; else, nu = max(nu, est); end
    slope = base - nu * rp1;
    sg = false;
    if ~(slope < 0) && rp1 > 0 && isfinite(base)
        nu = (base + 0.1 * abs(base) + realmin) / rp1; slope = base - nu * rp1; sg = true;
        F.nu_safeguard = F.nu_safeguard + 1;
    end
    last_dx = dx_inf;

    it = it_new(frame, k, F.main);
    it.tau_hat = tau / (Fw * L); it.rows = numel(rv); it.rd_inf = rd_inf / Fw; it.rp_inf_dhat = rp_inf / L;
    it.cent_hat = cent / (Fw * L); it.dx_tol = dx_inf / tol; it.min_gap_dhat = min(y) / dhat; it.E = E;
    it.inserted = inserted; it.levels = levels; it.nu_hat = nu / Fw; it.nu_safeguard = sg;
    if ~isempty(rv)
        slt = s .* lam / tau;
        it.slam_min = min(slt); it.slam_med = median(slt);
        cl = -Jdx ./ g; [it.closure_max, jc] = max(cl); it.closure_vertex = rv(jc);
    end
    if capture
        it.cap = struct('x', x, 'xt', xt, 'rv', rv, 's', s, 'lam', lam, 'src', src, 'tau', tau, 'g', g, 'dx', dx, ...
            'ds', ds, 'dlam', dlam, 'nu', nu, 'slope', slope, 'levels', levels, 'E', E, 'xk', xk);
    end

    % ---------------- stop test, outer event, final checks
    rconv = residuals_converged(tau, rp_inf, cent, dx_inf);
    conv = false;
    if rconv
        oc = so_outer('inner stop test');
        if oc == 2
            it.event = strjoin([ev, {'outer event failed'}], '; '); its = [its, it]; break; %#ok<AGROW>
        end
        if oc == 1
            F.discarded = F.discarded + 1; last_dx = -1;
            it.event = strjoin([ev, {'stop test -> outer event changed the state (system discarded)'}], '; ');
            it.limiter = 'discarded'; it.alpha = 0; it.alpha_eff = 0;
            its = [its, it]; %#ok<AGROW>
            continue;
        end
        full_ground = all(y + dy > 0);
        [first_full, np0] = elem_check(S, x, dx, 1);
        conv = full_ground && np0 == 0 && isinf(first_full);
        it.full_ccd = double(full_ground); it.full_valid = np0 == 0 && isinf(first_full);
        if ~conv, ev{end+1} = 'stop test passed, full-correction check failed'; end %#ok<AGROW>
    end
    terminal_now = Q.terminal_update && conv;
    if terminal_now
        F.terminal_candidates = F.terminal_candidates + 1; ev{end+1} = 'terminal candidate'; %#ok<AGROW>
    end

    % ---------------- step length: fractions to the boundary, update B / split rule, element validity
    [a_s, js] = ftb_min(s, ds, Q.ftb_eta);
    [a_l, jl] = ftb_min(lam, dlam, Q.ftb_eta);
    a_s = min(1, a_s); a_l = min(1, a_l);
    alpha = a_s; dual_lim = false;
    if a_l < alpha
        alpha = a_l; dual_lim = true; F.dual_limited = F.dual_limited + 1;
    end
    split_now = Q.split_tau1 && levels == 0 && a_l < a_s;
    if split_now
        alpha = a_s; F.split = F.split + 1;
    end
    [first, np] = elem_check(S, x, dx, alpha);
    if np > 0
        F.status = 'failed'; F.reason = 'an element is not valid at the start of the step';
        it.event = strjoin([ev, {F.reason}], '; '); its = [its, it]; break; %#ok<AGROW>
    end
    a_elem = Inf; lim_tri = 0;
    if first <= alpha
        a_elem = P.valid_eta * first; alpha = a_elem; F.valid_limited = F.valid_limited + 1; lim_tri = 1;
    end
    lim_hist{end+1} = lim_tri; %#ok<AGROW>
    a_prop = alpha;
    phi0 = E - tau * sum(log(s / dhat)) + nu * rp1;
    % ---------------- merit line search
    accepted = false; cert = 0; tr = 0;
    for tr = 1:P.ls_max
        xtr = x + alpha * dx;
        [~, npt] = elem_check(S, xtr);
        if npt > 0
            cert = cert + 1; alpha = alpha / 2; continue;
        end
        Etr = mech_eval(S, xtr, xt, 1);
        gtr = xtr(2*rv) - gy;
        str = max(s + alpha * ds, gtr);
        if split_now, a_lam = min(alpha, a_l); else, a_lam = alpha; end
        ltr = lam + a_lam * dlam;
        phi = Etr - tau * sum(log(str / dhat)) + nu * sum(abs(gtr - str));
        if armijo(phi0, phi, alpha, slope)
            accepted = true; break;
        end
        alpha = alpha / 2;
    end
    F.trials = F.trials + tr; F.cert_rejects = F.cert_rejects + cert; F.backtracks = F.backtracks + tr - 1;
    it.alpha_slack = a_s; it.alpha_mult = a_l; it.alpha_elem = a_elem; it.alpha_proposal = a_prop;
    it.split = split_now; it.trials = tr;
    if ~isempty(rv) && isfinite(a_s) && a_s < 1, it.lim_vertex = rv(js); end
    if ~accepted
        it.alpha = 0; it.alpha_eff = 0; it.limiter = 'line search rejected every trial';
        F.discarded = F.discarded + 1;
        if terminal_now
            F.terminal_rejected = F.terminal_rejected + 1;
            v = verify_now();
            it.event = strjoin([ev, {'terminal candidate rejected; verification of the unchanged state'}], '; ');
            its = [its, it, ver_trace(v)]; %#ok<AGROW>
            if v.hit
                F.status = 'converged'; F.reason = 'terminal (rejected) state verified';
            else
                F.status = 'failed'; F.reason = 'terminal update: candidate rejected and the unchanged state did not verify';
            end
            break;
        end
        F.status = 'failed';
        if cert == tr, F.reason = 'no trial position with valid elements'; else, F.reason = 'line search rejected every trial'; end
        it.event = strjoin([ev, {F.reason}], '; '); its = [its, it]; break; %#ok<AGROW>
    end
    % accept (positions, slacks, multipliers together)
    if split_now, a_lam = min(alpha, a_l); else, a_lam = alpha; end
    x = xtr; s = str; lam = ltr; x_moved = true;
    F.accepted = F.accepted + 1;
    it.alpha = alpha; it.alpha_lambda = a_lam; it.alpha_eff = alpha;
    if alpha < a_prop
        if cert == tr - 1 && cert > 0, it.limiter = 'element certificate'; else, it.limiter = 'merit backtracking'; end
    elseif isfinite(a_elem) && a_prop == a_elem
        it.limiter = 'element validity';
    elseif a_prop < 1
        if dual_lim && ~split_now
            it.limiter = 'multiplier positivity'; if ~isempty(rv), it.lim_vertex = rv(jl); end
        else
            it.limiter = 'slack positivity';
        end
    else
        it.limiter = 'none';
    end
    % the redirect trigger (two consecutive accepted steps cut to < 0.1 by an element): counted, not ported
    if numel(lim_hist) >= 2 && lim_hist{end} == 1 && lim_hist{end-1} == 1 && alpha < 0.1
        F.redirect_would_arm = F.redirect_would_arm + 1;
    end
    step = alpha * dx_inf;
    if step < Q.stall_ratio * tol, stalled = stalled + 1; else, stalled = 0; end
    it.event = strjoin(ev, '; ');
    its = [its, it]; %#ok<AGROW>
    if terminal_now
        F.terminal_applied = F.terminal_applied + 1;
        v = verify_now();
        its = [its, ver_trace(v)]; %#ok<AGROW>
        if v.hit
            F.status = 'converged'; F.reason = 'terminal update applied and verified';
            break;
        end
        F.terminal_continue = F.terminal_continue + 1;
        xk = x; so_refresh = true;   % a collision-safe outer position; refresh at the next iteration top
    end
    if stalled >= Q.stall_count
        F.status = 'failed'; F.reason = 'tiny steps'; break;
    end
end
if strcmp(F.status, 'running')
    F.status = 'budget'; F.reason = 'main-system budget reached';
end
F.levels = levels; F.rows_end = numel(rv); F.ver = ver_last;
if strcmp(F.status, 'converged')
    hist_out = struct('frame', frame, 'rv', rv, 'lam', lam);
end
st = struct('rv', rv, 's', s, 'lam', lam, 'src', src, 'tau', tau);

% ======================================================================================== nested helpers
    function v = JT(r)
        v = zeros(S.ndof, 1);
        v(2*rv) = r;
    end
    function tf = at_final(t)
        tf = t <= S.tau_final * (1 + 1e-12);
    end
    function tf = should_decrease(t, rdi, rpi, ce, ldx)
        tf = false;
        if at_final(t), return; end
        thh = t / (Fw * L);
        Et = max([rdi / Fw, rpi / L, ce / (Fw * L)]);
        if ~isfinite(Et), return; end
        if Et <= Q.kappa_eps * thh, tf = true; return; end
        tf = ldx >= 0 && ldx <= tol && rpi / L <= Q.kappa_eps * thh;
    end
    function tf = residuals_converged(t, rpi, ce, dxi)
        tf = at_final(t) && isfinite(dxi) && dxi <= tol && rpi <= tol && ce <= S.cent_final;
    end
    function tf = armijo(p0, p1, a, sl)
        tf = false;
        if ~isfinite(p1) || ~isfinite(p0), return; end
        if p1 <= p0 + Q.armijo_c * a * sl, tf = true; return; end
        tf = abs(Q.armijo_c * a * sl) < Q.precision * abs(p0) && p1 <= p0;
    end
    function add_rows(vids, code)
        if isempty(vids), return; end
        gg = x(2*vids) - gy;
        rv = [rv; vids(:)]; s = [s; gg]; lam = [lam; tau ./ gg];
        src = [src; code * ones(numel(vids), 1)]; born = [born; k * ones(numel(vids), 1)];
    end
    function [ok, nadd] = scan_or_fallback(trigger)
        % every vertex's signed ground gap at the accepted point
        ok = true; nadd = 0;
        yy = x(ydof) - gy;
        F.early_scans = F.early_scans + 1;
        if all(yy > 0)
            band = setdiff(find(yy < dhat), rv);
            if ~isempty(band)
                add_rows(band, 2);
                F.early_batches = F.early_batches + 1; F.early_rows = F.early_rows + numel(band);
                nadd = numel(band);
            end
            xk = x;
        else
            % a vertex at a nonpositive gap: the ground advance from the last checked position (no row is created)
            F.early_fallbacks = F.early_fallbacks + 1;
            F.fallback = F.fallback + 1;      % the system whose accepted step is clipped
            [oc, t_fb, nref] = outer_advance(['fallback (' trigger ')']);
            if oc == 2, ok = false; return; end
            nadd = nref;
            if ~isempty(its)
                j = find(strcmp({its.kind}, 'main'), 1, 'last');
                its(j).fallback_t = t_fb;
                its(j).alpha_eff = its(j).alpha * t_fb;
                its(j).limiter = [its(j).limiter ' + ground fallback'];
                its(j).event = strjoin({its(j).event, sprintf('ground fallback t=%.4g', t_fb)}, '; ');
            end
        end
    end
    function oc = so_outer(trigger)
        [oc, ~, ~] = outer_advance(trigger);
    end
    function [oc, t, nref] = outer_advance(trigger)
        % the outer ground advance from x_k toward x: full, then halved until every vertex gap is > 0; element bound
        oc = 0; t = 1; nref = 0;
        F.outer_events = F.outer_events + 1;
        seg = x - xk;
        yk = xk(ydof) - gy; sy = seg(ydof);
        tg = 1; h = 0;
        while ~all(yk + tg * sy > 0) && h < Q.ground_halvings
            tg = tg / 2; h = h + 1;
        end
        if ~all(yk + tg * sy > 0)
            F.status = 'failed'; F.reason = ['signed outer (ground): no strictly feasible fraction (' trigger ')'];
            oc = 2; return;
        end
        [fe, npk] = elem_check(S, xk, seg, 1);
        if npk > 0
            F.status = 'failed'; F.reason = 'signed outer: an element is not valid at the outer position'; oc = 2; return;
        end
        te = 1; if fe <= 1, te = P.valid_eta * fe; end
        t = min(tg, te);
        clipped = t < 1;
        if clipped
            okc = false;
            for hh = 1:Q.element_halvings
                xc = xk + t * seg;
                [~, npc] = elem_check(S, xc);
                if npc == 0, okc = true; break; end
                t = t / 2;
            end
            if ~okc
                F.status = 'failed'; F.reason = 'signed outer: no valid clipped position'; oc = 2; return;
            end
            x = xc;
            gg = x(2*rv) - gy;
            if any(gg <= 0)
                F.status = 'failed'; F.reason = 'signed outer: a row gap is not positive at the clipped position'; oc = 2; return;
            end
            s = gg;   % slacks set to the gaps of the motion actually taken; multipliers kept
            F.outer_clipped = F.outer_clipped + 1;
            adv = t * max(abs(seg));
            if adv < Q.stall_ratio * tol, consec_negl = consec_negl + 1; else, consec_negl = 0; end
            if consec_negl >= 3
                F.status = 'failed'; F.reason = 'signed outer: repeated clipping without progress'; oc = 2; return;
            end
        else
            consec_negl = 0;
        end
        xk = x;
        yy = x(ydof) - gy;
        band = setdiff(find(yy < dhat & yy > 0), rv);
        if ~isempty(band)
            add_rows(band, 3); F.refresh_rows = F.refresh_rows + numel(band); nref = numel(band);
        end
        x_moved = true; so_early_skip = true;
        if clipped || ~isempty(band), oc = 1; end
    end
    function v = verify_now()
        % the fresh verification of the returned state (one verification system)
        yy = x(ydof) - gy; gg = yy(rv);
        [~, gEv, H0v] = mech_eval(S, x, xt, 3);
        rdv = gEv - JT(lam);
        wv = lam ./ s;
        Kv = H0v + sparse(2*rv, 2*rv, wv, S.ndof, S.ndof);
        bv = gEv - JT(tau ./ s) + JT(wv .* (gg - s));
        dxv = shared_solve(Kv, -bv);
        F.verification = F.verification + 1;
        v = struct();
        v.rd_inf = max(abs(rdv)) / Fw;
        v.dx_tol = max(abs(dxv)) / tol;
        v.slack_consistency_dhat = maxabs(gg - s) / L;
        v.cent_hat = maxabs(s .* lam - tau) / (Fw * L);
        v.comp_dev = maxabs(s .* lam / tau - 1);
        v.at_final = at_final(tau);
        v.min_s = minval(s); v.min_lam = minval(lam);
        v.unresolved = numel(setdiff(find(yy < dhat), rv));
        v.ground_ok = all(yy > 0);
        v.min_gap_dhat = min(yy) / dhat;
        v.full_ground = all(yy + dxv(ydof) > 0);
        [fv, npv, dets] = elem_check(S, x, dxv, 1);
        v.elem_ok = npv == 0; v.full_valid = isinf(fv); v.min_det_rel = min(dets ./ (2 * S.A(:)));
        v.ok = v.dx_tol <= 1 && v.slack_consistency_dhat * L <= tol && v.cent_hat * Fw * L <= S.cent_final ...
            && v.at_final && v.min_s > 0 && v.min_lam > 0 && v.unresolved == 0 && v.ground_ok && v.full_ground ...
            && v.elem_ok && v.full_valid;
        v.hit = v.ok && v.rd_inf <= P.residual_target;
        v.solve = F.main;   % verification follows this many main systems
        ver_last = v;
    end
    function itv = ver_trace(v)
        itv = it_new(frame, k, F.main);
        itv.kind = 'verification'; itv.tau_hat = tau / (Fw * L); itv.rows = numel(rv);
        itv.rd_inf = v.rd_inf; itv.dx_tol = v.dx_tol; itv.rp_inf_dhat = v.slack_consistency_dhat;
        itv.cent_hat = v.cent_hat; itv.min_gap_dhat = v.min_gap_dhat;
        if v.hit, itv.event = 'verification passed'; else, itv.event = 'verification failed'; end
    end
end

% ======================================================================================== local functions
function it = it_new(frame, k, solve)
it = struct('solver', 'ip', 'frame', frame, 'k', k, 'solve', solve, 'kind', 'main', 'tau_hat', NaN, 'rows', 0, ...
    'rd_inf', NaN, 'rp_inf_dhat', NaN, 'cent_hat', NaN, 'dx_tol', NaN, 'alpha_slack', NaN, 'alpha_mult', NaN, ...
    'alpha_ccd', NaN, 'alpha_elem', NaN, 'alpha_proposal', NaN, 'alpha', NaN, 'alpha_lambda', NaN, 'fallback_t', NaN, ...
    'alpha_eff', NaN, 'trials', 0, 'limiter', '', 'lim_vertex', NaN, 'ccd_vertex', NaN, 'split', false, ...
    'native_stop', false, 'full_ccd', NaN, 'full_valid', NaN, 'min_gap_dhat', NaN, 'E', NaN, 'E_trial', NaN, ...
    'event', '', 'inserted', 0, 'levels', 0, 'nu_hat', NaN, 'nu_safeguard', false, 'slam_min', NaN, 'slam_med', NaN, ...
    'closure_max', NaN, 'closure_vertex', NaN, 'cap', []);
end

function [a, j] = ftb_min(v, dv, eta)
a = Inf; j = 0;
idx = find(dv < 0);
if isempty(idx), return; end
t = -eta * v(idx) ./ dv(idx);
[a, q] = min(t); j = idx(q);
end

function m = maxabs(v)
if isempty(v), m = 0; else, m = max(abs(v)); end
end

function m = minval(v)
if isempty(v), m = Inf; else, m = min(v); end
end
