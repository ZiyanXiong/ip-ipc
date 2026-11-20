function [name, row] = ip_step_limiter(SP, LS)
%IP_STEP_LIMITER  Which restriction set hard IP's accepted position step (the solver's own classification).
%   row: the contact row whose bound set the step (slack or multiplier positivity), 0 otherwise.
% The all-vertex ground safeguard acts after the accept (early-scan fallback) and is reported by the caller.
row = 0;
if ~LS.accepted
    name = 'line search rejected every trial';
elseif LS.alpha < SP.alpha_proposal
    if LS.cert == LS.trials - 1 && LS.cert > 0
        name = 'element certificate';
    else
        name = 'merit backtracking';
    end
elseif isfinite(SP.a_elem) && SP.alpha_proposal == SP.a_elem
    name = 'element validity';
elseif SP.alpha_proposal < 1
    if SP.dual_lim && ~SP.split_now
        name = 'multiplier positivity'; row = SP.jl;
    else
        name = 'slack positivity'; row = SP.js;
    end
else
    name = 'none';
end
end
