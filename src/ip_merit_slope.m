function [nu, slope, safeguard] = ip_merit_slope(base, rp1, lam_plus, nu_prev, Q)
%IP_MERIT_SLOPE  The l1 merit's weight and slope (ip_logic.h Params::merit_slope), shared by solver and diagnostics.
%   base = grad F . dx - tau sum(ds ./ s), rp1 = |g - s|_1, lam_plus = max|lambda + dlambda|.
% nu = nu_rho * lam_plus (a non-finite estimate falls back to max(nu_prev, estimate)); if the slope is not negative
% while |r_p|_1 > 0, the descent safeguard raises nu so that the slope is -(0.1 |base| + tiny).
est = Q.nu_rho * lam_plus;
if isfinite(est), nu = est; else, nu = max(nu_prev, est); end
slope = base - nu * rp1;
safeguard = false;
if ~(slope < 0) && rp1 > 0 && isfinite(base)
    nu = (base + 0.1 * abs(base) + realmin) / rp1;
    slope = base - nu * rp1;
    safeguard = true;
end
end
