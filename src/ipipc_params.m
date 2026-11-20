function P = ipipc_params()
%IPIPC_PARAMS  Every scene and solver setting of the 2D IP-vs-IPC box drop.
%
% Reference (read-only): E:\ZiyanXiong\github\libuipc-ip
%   scene      benchmarks/converged/eval_scenes.py  fem_block_drop (block 0.4 m, E 1e5, nu 0.4, rho 1e3,
%              clearance 0.3 m, tilt, frictionless ground at y = 0), dt 0.01, 200 frames
%   defaults   src/core/core/scene_default_config.cpp (d_hat 0.01, velocity_tol 0.05 m/s, line_search/max_iter 8,
%              newton/ccd_tol 1, cfl/enable 0, contact/adaptive/kappa_eval_scale 1e-16)
%   IP logic   src/backends/cuda/engine/ip_logic.h (Params defaults below)
%   IPC        src/backends/cuda/contact_system/global_contact_manager.cu (kappa corridor),
%              collision_detection/filters/easy_vertex_half_plane_trajectory_filter.cu (eta 0.1),
%              contact_models/sym/codim_ipc_contact.inl (KappaBarrierLog2)
%
% 2D normalization: the box is a plane-strain slab of unit depth (1 m). Masses, energies and forces are per metre
% of depth. The solver energy is E(x) = 1/2 |x - xt|_M^2 + dt^2 Psi_el(x) (+ IPC barrier), so its gradient has the
% unit force*dt^2 (kg m per m of depth). F_w = mean vertex mass * |g| * dt^2 is one vertex's weight in that unit.

P = struct();
P.name = 'ipipc 2D box drop (IP vs IPC), frictionless';

% ---- scene (2D analogue of fem_block_drop)
P.box_W      = 0.4;     % m
P.box_H      = 0.4;     % m
P.nx         = 16;      % cells along the width  (17 bottom vertices)
P.ny         = 16;      % cells along the height
P.tilt_deg   = 10;      % counter-clockwise; the bottom-left corner lands first
P.clearance  = 0.3;     % m, lowest vertex above the ground at t = 0
P.E          = 1e5;     % Pa
P.nu         = 0.4;
P.rho        = 1e3;     % kg/m^3 (x 1 m depth -> kg/m^2 areal density)
P.gravity    = [0; -9.8];
P.dt         = 0.01;    % s
P.frames     = 200;
P.ground_y   = 0.0;     % static line y = 0, normal (0, 1)
P.material   = 'stable neo-Hookean 2D (rest-stable form, plane-strain Lame parameters)';

% ---- shared numerics
P.d_hat           = 0.01;   % m: IPC activation distance and IP discovery band (contact/d_hat default)
P.velocity_tol    = 0.05;   % m/s (newton/velocity_tol default); tol_disp = velocity_tol * dt
P.residual_target = 5;      % ||r_d||_inf target in F_w at the returned state (ip/residual_target default)
P.ls_max          = 8;      % line_search/max_iter (both solvers)
P.cap_main        = 256;    % main Newton systems per frame (common budget)
P.valid_eta       = 0.9;    % element validity: keep 90% of the first inversion time (ip_logic.h valid_eta)

% ---- IPC (native ground contact)
P.ipc.kappa_default    = 1e9;     % fem_block_drop: contact_tabular().default_model(mu, 1e9)
P.ipc.kappa_eval_scale = 1e-16;   % contact/adaptive/kappa_eval_scale
P.ipc.min_kappa_coef   = 1e11;    % Stiff-GIPC minKappaCoef
P.ipc.corridor_factor  = 100;     % upper bound = 100 x lower bound
P.ipc.ccd_eta          = 0.1;     % half-plane CCD keeps (1 - eta) of the gap
P.ipc.ccd_tol          = 1.0;     % newton/ccd_tol: full-correction CCD fraction required at the matched stop
P.ipc.cfl              = false;   % cfl/enable 0

% ---- hard IP, round-43 'base' arm (arm H + signed schedule + zero-CCD ground advance + early discovery)
P.ip.tau0_hat        = 1.0;    % tau0 = F_w L
P.ip.tau_final_hat   = -1;     % <= 0: tol_disp / L
P.ip.kappa_eps       = 10.0;   % subproblem tolerance E_tau <= kappa_eps tau_hat; centrality bound kappa_eps tau_final
P.ip.kappa_tau       = 0.2;    % next_tau (unused after the direct first decrease, kept for completeness)
P.ip.theta_tau       = 1.5;
P.ip.ftb_eta         = 0.995;  % slack and multiplier fraction to the boundary
P.ip.s_floor_hat     = 1e-8;   % initial slack floor, in L
P.ip.armijo_c        = 1e-4;
P.ip.nu_rho          = 1.1;    % nu = nu_rho max|lambda + dlambda|
P.ip.precision       = 1e-12;
P.ip.stall_ratio     = 1e-6;
P.ip.stall_count     = 3;
P.ip.history         = true;   % UIPC_IP_MULTIPLIER_HISTORY=1
P.ip.terminal_update = true;   % UIPC_IP_TERMINAL_UPDATE=1
P.ip.split_tau1      = true;   % split multiplier rule at tau_hat = 1
P.ip.tau_direct      = true;   % first decrease goes directly to the final level
P.ip.early_discovery = true;   % UIPC_IP_SIGNED_EARLY_DISCOVERY=1
P.ip.ground_advance  = true;   % UIPC_IP_SIGNED_OUTER=1, _ADVANCE=ground (no CCD; halving to strictly positive gaps)
P.ip.ground_halvings = 40;
P.ip.element_halvings = 20;
P.ip.not_ported = {'redirected step (element-collapse trigger; never needed if no element limits two steps in a row: counted)', ...
                   'reduced first step, P70 buffer, logarithmic gap update (left for future experiments)', ...
                   'self contact (not modeled in 2D)'};

% ---- impact-frame selection for the replay
P.short_step_threshold = 0.25;  % "clearly shortened": accepted position fraction <= 0.25 of the Newton step
P.impact_frames        = 20;    % impact window length after first contact (reference: frames 24-43)
end
