function SP = ip_step_proposal(S, P, x, s, lam, ds, dlam, dx, levels)
%IP_STEP_PROPOSAL  Hard IP's proposed position fraction before the line search, shared by solver and diagnostics.
% Slack and multiplier fractions to the boundary (eta); update B (one common fraction) except at the first barrier
% level (levels == 0), where the split rule proposes the slack fraction for positions and slacks; then the element
% bound (valid_eta of the first inversion time on (0, alpha]).
%   SP.a_s / SP.js   slack fraction (<= 1) and its row (0: no slack decreases)
%   SP.a_l / SP.jl   multiplier fraction (<= 1) and its row (0: no multiplier decreases)
%   SP.np            elements not valid at the start (the solver fails the frame if > 0)
%   SP.a_elem, SP.elem_tri   the element bound (Inf if none) and the triangle that sets it
Q = P.ip;
[a_s, js] = ftb_min(s, ds, Q.ftb_eta);
[a_l, jl] = ftb_min(lam, dlam, Q.ftb_eta);
a_s = min(1, a_s); a_l = min(1, a_l);
alpha = a_s; dual_lim = false;
if a_l < alpha
    alpha = a_l; dual_lim = true;
end
split_now = Q.split_tau1 && levels == 0 && a_l < a_s;
if split_now
    alpha = a_s;
end
alpha_ftb = alpha;
[first, np, ~, tri] = elem_check(S, x, dx, alpha);
a_elem = Inf; elem_lim = false; elem_tri = NaN;
if np == 0 && first <= alpha
    a_elem = P.valid_eta * first; alpha = a_elem; elem_lim = true; elem_tri = tri;
end
SP = struct('a_s', a_s, 'js', js, 'a_l', a_l, 'jl', jl, 'dual_lim', dual_lim, 'split_now', split_now, ...
    'alpha_ftb', alpha_ftb, 'np', np, 'first', first, 'a_elem', a_elem, 'elem_lim', elem_lim, 'elem_tri', elem_tri, ...
    'alpha_proposal', alpha);
end

function [a, j] = ftb_min(v, dv, eta)
a = Inf; j = 0;
idx = find(dv < 0);
if isempty(idx), return; end
t = -eta * v(idx) ./ dv(idx);
[a, q] = min(t); j = idx(q);
end
