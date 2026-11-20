function [ok, by_precision] = ip_armijo(p0, p1, a, sl, Q)
%IP_ARMIJO  Hard IP's line-search acceptance (ip_logic.h Params::armijo), shared by the solver and the diagnostics.
% Armijo at the actual trial; when the predicted decrease is below the energy's rounding (precision * |phi0|), a
% non-increase is accepted (the roundoff rule); a non-finite merit is never accepted.
ok = false; by_precision = false;
if ~isfinite(p1) || ~isfinite(p0), return; end
if p1 <= p0 + Q.armijo_c * a * sl
    ok = true; return;
end
if abs(Q.armijo_c * a * sl) < Q.precision * abs(p0) && p1 <= p0
    ok = true; by_precision = true;
end
end
