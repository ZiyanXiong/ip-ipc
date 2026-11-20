function R = check_derivatives(S, P)
%CHECK_DERIVATIVES  Numerical checks of the IPC barrier's first and second derivatives at representative positive
% gaps (scalar and through the vertex-position assembly), and of the shared mechanical gradient and Hessian.
Dhat = S.d_hat^2; k = S.kappa_dt2;
gaps = [0.02 0.05 0.1 0.25 0.5 0.75 0.9 0.99]' * S.d_hat;
rows = [];
for q = 1:numel(gaps)
    d = gaps(q);
    % (a) in D: B'(D), B''(D) against central differences of B and B'
    hD = 1e-6 * d^2;
    [B0, dB0, ddB0] = ipc_barrier(d^2, Dhat); %#ok<ASGLU>
    [Bp, dBp] = ipc_barrier(d^2 + hD, Dhat); [Bm, dBm] = ipc_barrier(d^2 - hD, Dhat);
    e_dB = abs((Bp - Bm) / (2*hD) - dB0) / abs(dB0);
    e_ddB = abs((dBp - dBm) / (2*hD) - ddB0) / abs(ddB0);
    % (b) through the vertex position: energy -> gradient -> Hessian (the solver's own assembly)
    x = S.x0; x(2:2:end) = x(2:2:end) + 1;            % lift everything out of the band
    v = S.corner_first; x(2*v) = S.ground_y + d;
    [E0, g0, H0] = ipc_contact_eval(S, x, 3);
    h = 1e-6 * d;
    xp = x; xp(2*v) = xp(2*v) + h; xm = x; xm(2*v) = xm(2*v) - h;
    [Ep, gp] = ipc_contact_eval(S, xp, 3); [Em, gm] = ipc_contact_eval(S, xm, 3);
    fd_g = (Ep - Em) / (2*h); an_g = g0(2*v);
    fd_H = (gp(2*v) - gm(2*v)) / (2*h); an_H = H0(2*v, 2*v);
    h2 = 1e-3 * d;
    xp2 = x; xp2(2*v) = xp2(2*v) + h2; xm2 = x; xm2(2*v) = xm2(2*v) - h2;
    fd_H2 = (ipc_contact_eval(S, xp2, 1) - 2*E0 + ipc_contact_eval(S, xm2, 1)) / h2^2;
    rows = [rows; d / S.d_hat, -an_g, an_H, abs(fd_g - an_g) / abs(an_g), abs(fd_H - an_H) / abs(an_H), ...
        abs(fd_H2 - an_H) / abs(an_H), e_dB, e_ddB, -an_g / S.F_w, an_H / S.F_w * S.d_hat]; %#ok<AGROW>
end
R.barrier = array2table(rows, 'VariableNames', {'gap_over_dhat', 'force_fn', 'stiffness_kn', 'relerr_grad_fd', ...
    'relerr_hess_fd_of_grad', 'relerr_hess_2nd_fd_energy', 'relerr_dBdD', 'relerr_d2BdD2', 'force_over_Fw', ...
    'stiffness_times_dhat_over_Fw'});
% the stiffness stays positive over the band (the half-plane Hessian is not projected natively)
dd = linspace(1e-4, 0.9999, 2000)' * S.d_hat;
[~, dB, ddB] = ipc_barrier(dd.^2, Dhat);
kn = k * (4 * dd.^2 .* ddB + 2 * dB);
R.barrier_min_kn_over_band = min(kn);
R.barrier_kn_positive = all(kn > 0);

% shared mechanics: gradient by directional differences, the unprojected Hessian by differences of the gradient
rng(7);
hel = P.box_W / P.nx;
x = S.x0 + 0.1 * hel * randn(S.ndof, 1);
xt = S.x0 + 0.05 * hel * randn(S.ndof, 1);
[E0, g0, H0] = mech_eval(S, x, xt, 3, false);
p = randn(S.ndof, 1);
h = 1e-6;
Ep = mech_eval(S, x + h*p, xt, 1); Em = mech_eval(S, x - h*p, xt, 1);
[~, gp] = mech_eval(S, x + h*p, xt, 2); [~, gm] = mech_eval(S, x - h*p, xt, 2);
R.mech_grad_relerr = abs((Ep - Em) / (2*h) - g0' * p) / abs(g0' * p);
R.mech_hess_relerr = norm((gp - gm) / (2*h) - H0 * p) / norm(H0 * p);
[~, ~, Hp] = mech_eval(S, x, xt, 3, true);
R.mech_projected_min_eig = eigs(Hp, 1, 'smallestreal');
R.mech_energy_at_rest = mech_eval(S, S.x0, S.x0, 1);
end
