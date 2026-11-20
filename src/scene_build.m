function S = scene_build(P)
%SCENE_BUILD  Shared 2D scene: triangulated tilted box, lumped masses, material, units, IPC kappa.
nx = P.nx; ny = P.ny; W = P.box_W; H = P.box_H;
[I, J] = meshgrid(0:nx, 0:ny);            % I along x, J along y
I = I'; J = J';                            % vertex id = j*(nx+1) + i + 1
X0 = [-W/2 + I(:)' * W / nx; -H/2 + J(:)' * H / ny];
vid = @(i, j) j * (nx + 1) + i + 1;
T = zeros(2 * nx * ny, 3); t = 0;
for j = 0:ny-1
    for i = 0:nx-1
        a = vid(i, j); b = vid(i+1, j); c = vid(i+1, j+1); d = vid(i, j+1);
        if mod(i + j, 2) == 0
            T(t+1, :) = [a b c]; T(t+2, :) = [a c d];
        else
            T(t+1, :) = [a b d]; T(t+2, :) = [b c d];
        end
        t = t + 2;
    end
end
th = deg2rad(P.tilt_deg);
R = [cos(th) -sin(th); sin(th) cos(th)];
X = R * X0;
X(2, :) = X(2, :) + P.clearance - min(X(2, :));
X(1, :) = X(1, :) - mean(X(1, :));
N = size(X, 2);

% rest shape = initial (rigidly rotated) shape
x1 = X(:, T(:,1)); x2 = X(:, T(:,2)); x3 = X(:, T(:,3));
Dm11 = x2(1,:) - x1(1,:); Dm21 = x2(2,:) - x1(2,:);
Dm12 = x3(1,:) - x1(1,:); Dm22 = x3(2,:) - x1(2,:);
detDm = Dm11 .* Dm22 - Dm12 .* Dm21;
assert(all(detDm > 0), 'triangles must be counter-clockwise');
A = 0.5 * detDm;
B11 =  Dm22 ./ detDm; B12 = -Dm12 ./ detDm;
B21 = -Dm21 ./ detDm; B22 =  Dm11 ./ detDm;   % B = inv(Dm)

m = accumarray(T(:), repmat(P.rho * A(:) / 3, 3, 1), [N 1]);   % lumped, kg per m of depth

S = struct();
S.N = N; S.ndof = 2 * N; S.T = T; S.nT = size(T, 1);
S.X0 = X; S.x0 = X(:);
S.A = A; S.B11 = B11; S.B12 = B12; S.B21 = B21; S.B22 = B22;
S.m = m; S.m2 = kron(m, [1; 1]);
S.mu  = P.E / (2 * (1 + P.nu));
S.lam = P.E * P.nu / ((1 + P.nu) * (1 - 2 * P.nu));
S.dt = P.dt; S.dt2 = P.dt^2;
S.gvec = repmat(P.gravity, N, 1);
S.ground_y = P.ground_y;
S.ydof = (2:2:2*N)';
S.dofs = [2*T(:,1)-1, 2*T(:,1), 2*T(:,2)-1, 2*T(:,2), 2*T(:,3)-1, 2*T(:,3)];   % nT x 6
% element vec(F) Jacobian (vec order F11 F21 F12 F22; xe = [x1 y1 x2 y2 x3 y3]), 4 x 6 x nT
D = zeros(4, 6, S.nT);
cx1 = -(B11 + B21); cx2 = B11; cx3 = B21;   % F_i1 coefficients
cy1 = -(B12 + B22); cy2 = B12; cy3 = B22;   % F_i2 coefficients
D(1, 1, :) = cx1; D(1, 3, :) = cx2; D(1, 5, :) = cx3;   % F11 (x comps)
D(2, 2, :) = cx1; D(2, 4, :) = cx2; D(2, 6, :) = cx3;   % F21 (y comps)
D(3, 1, :) = cy1; D(3, 3, :) = cy2; D(3, 5, :) = cy3;   % F12
D(4, 2, :) = cy1; D(4, 4, :) = cy2; D(4, 6, :) = cy3;   % F22
S.Dmat = D;
[r, c] = ndgrid(1:6, 1:6);
S.Hrow = reshape(S.dofs(:, r(:))', [], 1);
S.Hcol = reshape(S.dofs(:, c(:))', [], 1);

% bottom-edge vertices (reference configuration)
S.bottom = find(J(:) == 0);
S.corner_first = vid(0, 0);

% ---- units (see ipipc_params.m)
S.mbar = mean(m);
S.F_w  = S.mbar * norm(P.gravity) * P.dt^2;     % solver force unit
S.L    = P.d_hat;
S.d_hat = P.d_hat;
S.tol_disp = P.velocity_tol * P.dt;             % m
S.tau0 = P.ip.tau0_hat * S.F_w * S.L;
if P.ip.tau_final_hat > 0
    S.tau_final = P.ip.tau_final_hat * S.F_w * S.L;
else
    S.tau_final = S.tol_disp / S.L * S.F_w * S.L;
end
S.cent_final = P.ip.kappa_eps * S.tau_final;
S.s_floor = P.ip.s_floor_hat * S.L;

% ---- IPC kappa: the native scene-adaptive corridor (global_contact_manager.cu, Stiff-GIPC port), in 2D.
% The scene diagonal is the vertex bounding box of the rest scene, including the ground's reference point.
allX = [X, [0; P.ground_y]];
diagS = norm(max(allX, [], 2) - min(allX, [], 2));
mean_mass = sum(m) / (N + 1);                      % total mass / total vertex count (ground vertex included)
d = P.ipc.kappa_eval_scale * diagS^2;
dHat = P.d_hat^2;
tt = d - dHat;
H_b = -2 * log(d / dHat) - 4 * tt / d + tt^2 / d^2;
base = P.ipc.min_kappa_coef / (4 * d * H_b) * mean_mass;
k_lo = base / P.dt^2;
k_hi = P.ipc.corridor_factor * k_lo;
S.kappa = min(max(P.ipc.kappa_default, k_lo), k_hi);
S.kappa_corridor = [k_lo, k_hi];
S.kappa_dt2 = S.kappa * P.dt^2;      % the coefficient of the barrier in the solver energy
S.scene_diagonal = diagS;
S.kappa_mean_mass = mean_mass;
end
