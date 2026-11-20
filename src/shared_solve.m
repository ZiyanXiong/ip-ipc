function dx = shared_solve(K, rhs)
%SHARED_SOLVE  The common full-size sparse linear solver of both methods (sparse Cholesky, CHOLMOD via chol).
% Every call is one full-size solve; the callers count them.
[R, p, Q] = chol(K, 'vector');
if p ~= 0
    % not positive definite (should not happen: PSD-projected elasticity + mass + nonnegative contact terms)
    dx = K \ rhs;
    return;
end
dx = zeros(size(rhs));
dx(Q, :) = R \ (R' \ rhs(Q, :));
end
