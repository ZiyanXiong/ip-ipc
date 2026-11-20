function LS = ip_line_search(S, P, x, xt, rv, s, lam, dx, ds, dlam, tau, nu, slope, phi0, alpha, split_now, a_l)
%IP_LINE_SEARCH  Hard IP's merit line search, shared by the solver and the diagnostics.
% Trials x~ = x + alpha dx (halving, at most ls_max): every element must be valid at the trial (else halve); then
%   s~ = max(s + alpha ds, g(x~)),  lambda~ = lambda + alpha_lambda dlambda (alpha_lambda = min(alpha, a_l) under
%   the split rule, alpha otherwise),  phi = E_mech(x~) - tau sum log(s~/d_hat) + nu sum|g(x~) - s~|
% accepted by ip_armijo (Armijo, or the roundoff rule).
%   LS.log rows: [alpha, phi - phi0, armijo_c alpha slope, accepted, by_precision, element_rejected]
Q = P.ip;
LS = struct('accepted', false, 'alpha', alpha, 'a_lam', NaN, 'trials', 0, 'cert', 0, 'by_precision', false, ...
    'x', [], 's', [], 'lam', [], 'phi', NaN, 'log', zeros(0, 6));
cert = 0; tr = 0;
for tr = 1:P.ls_max
    xtr = x + alpha * dx;
    [~, npt] = elem_check(S, xtr);
    if npt > 0
        LS.log(end+1, :) = [alpha, NaN, NaN, 0, 0, 1];
        cert = cert + 1; alpha = alpha / 2; continue;
    end
    Etr = mech_eval(S, xtr, xt, 1);
    gtr = xtr(2*rv) - S.ground_y;
    str = max(s + alpha * ds, gtr);
    if split_now, a_lam = min(alpha, a_l); else, a_lam = alpha; end
    ltr = lam + a_lam * dlam;
    phi = Etr - tau * sum(log(str / S.d_hat)) + nu * sum(abs(gtr - str));
    [ok, byp] = ip_armijo(phi0, phi, alpha, slope, Q);
    LS.log(end+1, :) = [alpha, phi - phi0, Q.armijo_c * alpha * slope, ok, byp, 0];
    if ok
        LS.accepted = true; LS.by_precision = byp;
        LS.x = xtr; LS.s = str; LS.lam = ltr; LS.a_lam = a_lam; LS.phi = phi;
        break;
    end
    alpha = alpha / 2;
end
LS.alpha = alpha; LS.trials = tr; LS.cert = cert;
end
