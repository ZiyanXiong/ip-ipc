function [Eb, gb, Hb, C] = ipc_contact_eval(S, x, want)
%IPC_CONTACT_EVAL  Native IPC vertex-ground barrier at positions x (ipc_vertex_half_plane_contact_function.h).
% The active pairs are evaluated at x itself: every vertex with D = g^2 < d_hat^2 (and g > 0).
%   energy   kappa dt^2 B(D)
%   gradient kappa dt^2 B'(D) dD/dx,            dD/dx = 2 g n
%   Hessian  kappa dt^2 (B''(D) dD dD' + B'(D) d2D/dx2),  d2D/dx2 = 2 n n'   (not projected, as natively)
% Normal quantities in solver units: force f_n = -d(kappa dt^2 B)/dg = -kappa dt^2 B'(D) 2 g,
% stiffness k_n = d2(kappa dt^2 B)/dg2 = kappa dt^2 (4 g^2 B''(D) + 2 B'(D)).
if nargin < 3, want = 3; end
y = x(S.ydof) - S.ground_y;
D = y.^2;
Dhat = S.d_hat^2;
act = find(y < S.d_hat & y > 0);
C.act = act; C.g = y(act);
[B, dB, ddB] = ipc_barrier(D(act), Dhat);
k = S.kappa_dt2;
Eb = k * sum(B);
if any(y <= 0)
    Eb = Inf;  % a vertex at or below the ground: infeasible for IPC
end
C.fn = -k * dB .* 2 .* y(act);
C.kn = k * (4 * y(act).^2 .* ddB + 2 * dB);
gb = []; Hb = [];
if want < 2, return; end
gb = zeros(S.ndof, 1);
gb(2*act) = -C.fn;
if want < 3, return; end
Hb = sparse(2*act, 2*act, C.kn, S.ndof, S.ndof);
end
