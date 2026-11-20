function [E, g, H, aux] = mech_eval(S, x, xt, want, project)
%MECH_EVAL  Shared mechanical energy and derivatives (both solvers use exactly this).
%   E(x) = 1/2 (x - xt)' M (x - xt) + dt^2 sum_e A_e Psi(F_e)
%   Psi  = mu/2 (|F|^2 - 2) - mu (J - 1) + lam/2 (J - 1)^2     (rest-stable 2D stable neo-Hookean)
%   want: 1 energy, 2 + gradient, 3 + Hessian (per-element PSD projection of d2Psi/dF2 unless project == false)
if nargin < 4, want = 3; end
if nargin < 5, project = true; end
T = S.T;
X1x = x(2*T(:,1)-1)'; X1y = x(2*T(:,1))';
Ds11 = x(2*T(:,2)-1)' - X1x; Ds21 = x(2*T(:,2))' - X1y;
Ds12 = x(2*T(:,3)-1)' - X1x; Ds22 = x(2*T(:,3))' - X1y;
F11 = Ds11 .* S.B11 + Ds12 .* S.B21;  F12 = Ds11 .* S.B12 + Ds12 .* S.B22;
F21 = Ds21 .* S.B11 + Ds22 .* S.B21;  F22 = Ds21 .* S.B12 + Ds22 .* S.B22;
Jd = F11 .* F22 - F12 .* F21;
mu = S.mu; lam = S.lam;
Psi = mu/2 * (F11.^2 + F12.^2 + F21.^2 + F22.^2 - 2) - mu * (Jd - 1) + lam/2 * (Jd - 1).^2;
w = S.dt2 * S.A;
r = x - xt;
E = 0.5 * sum(S.m2 .* r.^2) + sum(w .* Psi);
aux.J = Jd;
g = []; H = [];
if want < 2, return; end
c = lam * (Jd - 1) - mu;
P11 = mu * F11 + c .* F22;   P12 = mu * F12 - c .* F21;
P21 = mu * F21 - c .* F12;   P22 = mu * F22 + c .* F11;
g2x = P11 .* S.B11 + P12 .* S.B12;  g2y = P21 .* S.B11 + P22 .* S.B12;
g3x = P11 .* S.B21 + P12 .* S.B22;  g3y = P21 .* S.B21 + P22 .* S.B22;
ge = w .* [-(g2x + g3x); -(g2y + g3y); g2x; g2y; g3x; g3y];      % 6 x nT
g = S.m2 .* r + accumarray(reshape(S.dofs', [], 1), ge(:), [S.ndof 1]);
if want < 3, return; end
nT = S.nT;
cv = reshape([F22; -F12; -F21; F11], 4, 1, nT);                  % vec(cof F) (vec order F11 F21 F12 F22)
HF = repmat(mu * eye(4), 1, 1, nT) + lam * pagemtimes(cv, 'none', cv, 'transpose');
cc = reshape(c, 1, 1, nT);
HF(1,4,:) = HF(1,4,:) + cc; HF(4,1,:) = HF(4,1,:) + cc;
HF(2,3,:) = HF(2,3,:) - cc; HF(3,2,:) = HF(3,2,:) - cc;
if project
    HF = 0.5 * (HF + pagetranspose(HF));
    [V, Dg] = pageeig(HF);
    Dg = max(real(Dg), 0);
    V = real(V);
    HF = pagemtimes(pagemtimes(V, Dg), 'none', V, 'transpose');
end
He = pagemtimes(pagemtimes(S.Dmat, 'transpose', HF, 'none'), S.Dmat);   % 6 x 6 x nT
He = He .* reshape(w, 1, 1, nT);
H = sparse(S.Hrow, S.Hcol, He(:), S.ndof, S.ndof) + spdiags(S.m2, 0, S.ndof, S.ndof);
H = 0.5 * (H + H');
end
