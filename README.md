# ipipc-matlab — 2D box drop: hard IP vs IPC Newton solves

A standalone MATLAB port of the `libuipc-ip` block-drop comparison to 2D. A tilted elastic box falls onto a frictionless ground line and is simulated twice: by native IPC (log² barrier, CCD, energy line search, matched stop) and by ordinary hard interior point (IP; round 43's `base` arm). Both share the mesh, masses, material, time step, inertial prediction, element checks and sparse solver. The question the project is set up for is why hard IP needs more Newton solves than IPC. Findings so far, and what is still open, are in [REPORT.md](REPORT.md).

## Run it

Tested with **MATLAB R2025b Update 1 (25.2.0.3042426)** on Windows 11. Only base MATLAB is required (no toolboxes). From the project folder:

```
matlab -batch "run_experiment"
```

This takes about 15 s plus MATLAB startup. Every invocation writes a **new** folder, `runs/run_<yyyyMMdd_HHmmss>_<id>/`; existing folders and anything under `reference/` are refused. To choose the folder, use `run_experiment('OutDir', 'path/that/does/not/exist')`.

```
matlab -batch "run_tests"                                               % focused checks, about 30 s
matlab -batch "addpath tools; compare_to_reference('runs/run_...')"     % a run against the reference
matlab -batch "addpath tools; verify_reference"                         % reference integrity (SHA-256)
```

## Expected results (identical on the tested version: the code is deterministic)

| quantity | hard IP | IPC |
|---|---|---|
| qualified frames | 200/200 | 200/200 |
| full run: main + verification = total full-size solves | 415 + 200 = 615 | 362 + 0 = 362 |
| impact, frames 24–43 | 80 + 20 = 100 | 67 |
| free flight 1–23 / settling 44–200 (main) | 46 / 289 | 46 / 249 |
| max ‖r_d‖∞ at returned states (target 5 F_w) | 0.122 F_w | 4.906 F_w |
| replays of the 20 impact frames reproducing their trajectory frame | 20/20 | 20/20 |

`compare_to_reference.m` checks all of these against the preserved run, frame by frame, and names the first frame and quantity that differ.

## Replay a saved frame

Every run saves each solver's frame-start states (positions, velocities, inertial prediction and, for IP, the multipliers carried into the frame) in `trajectories.mat`:

```matlab
addpath tools
[F, its] = replay_frame('reference/baseline_2026-11-20/results', 26, 'ip');   % or a runs/run_... folder; 'ipc'
its(1).cap        % first iterate: x, rows, s, lambda, tau, dx, ds, dlambda, nu, ...
[F, its, D] = replay_frame('reference/baseline_2026-11-20/results', 26, 'ip', true);   % + diagnose_impact
D.limiting        % first limiting IP iterate: actual limiter, restrictions, centered-multiplier and IPC directions
```

The same can be done by hand:

```matlab
addpath src
P = ipipc_params(); S = scene_build(P);
M = load('reference/baseline_2026-11-20/results/trajectories.mat', 'Rip');
st = M.Rip.starts{26};
[x, F, its] = ip_solve_frame(S, P, st.x, st.xt, 26, st.hist, true);
```

`diagnose_iterate(S, P, its(k))` analyses any captured iterate. It re-forms the step with the solver's own shared rules (`ip_step_proposal`, `ip_line_search`, `ip_armijo`, `ip_step_limiter`) and reports whichever restriction actually set the step.

## Layout

| path | role |
|---|---|
| `run_experiment.m` | entry point: both trajectories, impact-frame replays, outputs |
| `run_tests.m`, `tests/` | diagnostics on real iterates, the controlled-failure path, runner guards, the reference manifest |
| `src/ipipc_params.m` | every scene and solver setting |
| `src/ip_solve_frame.m`, `src/ip_*.m` | hard IP (one frame) and its step rules, shared with the diagnostics |
| `src/ipc_solve_frame.m`, `src/ipc_*.m` | IPC (one frame), the native barrier |
| `src/mech_eval.m`, `elem_check.m`, `shared_solve.m`, `scene_build.m` | shared mechanics, element validity, sparse solver, scene |
| `src/run_trajectory.m`, `compare_performance.m` | trajectories (stop at the first failed frame); comparison only of complete trajectories |
| `src/diagnose_impact.m`, `diagnose_iterate.m` | same-state replays, replay agreement, first-limiting-iterate analysis |
| `tools/` | `compare_to_reference`, `verify_reference`, `replay_frame` |
| `reference/` | the preserved verified run, its code, SHA-256 manifest, and the identity of the 3D sources ([reference/README.md](reference/README.md)) |
| `runs/` | new runs (git-ignored) |

**Failed frames.** A frame qualifies only if its solver's stop test and returned-state checks pass. The first frame that fails ends that trajectory:
- its diagnostics are saved (`failure_<solver>.json`, its iterations, its unqualified returned state);
- position, velocity, multiplier history and time stay at the last qualified frame;
- `performance_comparison.json` reports "not comparable" instead of comparing solve counts or wall time.

## The 3D reference

The 3D reference is the `fem_block_drop` study in `E:\ZiyanXiong\github\libuipc-ip` (read-only): round 20 (early ground discovery), round 21 (equations, step limits) and round 43's `base` arm (the hard-IP configuration). That checkout is at commit `b611ae72`, but the hard-IP solver files are untracked and several others carry uncommitted edits, so the commit does not identify them. `reference/libuipc_ip_sources.sha256` records the SHA-256 of every file the port was read from. Its `ip_solver.cu` digest equals the one in round 43's `identity.txt`. Details: [reference/README.md](reference/README.md).

## Moving a result from 2D to 3D: differences that matter

- **Contact topology.**
  - 2D: 289 vertices, 17 on the bottom edge. The box is tilted 10° about one axis, lands on a corner and rotates onto a 1D contact edge.
  - 3D: 343 vertices, 49 on the bottom face, tilted about two axes (15°, 10°). It lands on a vertex onto a 2D patch that can rock in two directions.
  - The coupled "near rows open, far rows close" motion behind the short IP steps has one rotation axis in 2D and two in 3D; row counts per frame differ (≤ 17 against ≤ 49).
- **Units.**
  - Force tolerances are in F_w (one vertex's weight × dt²), so the 5 F_w target and τ̂ = 1 → 0.05 mean the same per vertex.
  - The body weighs 289 F_w in 2D and 343 F_w in 3D, and the load per contact row differs.
  - 2D quantities are per metre of depth (plane strain).
- **Material.** The 2D model is a rest-stable stable-neo-Hookean form with plane-strain Lamé parameters, not libuipc's 3D StableNeoHookean. With the same E, ν and ρ, impact deformation differs.
- **IPC stiffness.** κ uses the native corridor rule, with the 2D scene diagonal (0.893 m) and mean mass. κ·dt²·d̂⁴/F_w is 0.203 in 2D against about 0.30 in 3D, so IPC's barrier is relatively softer here.
- **Linear algebra and time.** Exact sparse Cholesky stands in for GPU PCG. Solve counts transfer; PCG iterations and wall time do not, and MATLAB timings say nothing about GPU timings.
- **Not ported:**
  - IP's element-redirect step (its trigger never armed here);
  - self contact;
  - the rounding-certified determinant test (det > 0 is used instead);
  - friction;
  - the reduced first step, P70 buffer and logarithmic gap update.
- **Orientation.** In 3D round 43, `base` used 547 main + 200 verification against IPC's 530 (matched stop), and 107 + 20 against 94 in impact. Here: 415 + 200 against 362, and 80 + 20 against 67.
