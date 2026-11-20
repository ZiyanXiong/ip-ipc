# 2D box drop: hard IP vs IPC — does the carried-multiplier short step reproduce, and what does it cost?

Standalone MATLAB experiment. Reference (read only): `E:\ZiyanXiong\github\libuipc-ip` (round 20 for discovery, round 21 for the equations and step limits, round 43 `base` for the hard-IP configuration). All code and results are in this folder.

**Reproduce** (about 12 s plus MATLAB startup; R2025b):

```powershell
matlab -batch "cd('E:\ZiyanXiong\github\ipipc-matlab'); run_experiment"
```

The run is deterministic. Every replay reproduces its trajectory frame exactly (zero difference in step fractions and residuals).

## Answers

1. **Did the short-step mechanism reproduce? Yes, but stiffness is only part of the cause.** In every short-step impact frame (25–28) the first Newton iterate runs at τ̂ = 1, the multipliers are carried from the previous frame, and every one of them sits at sλ/μ ≈ 0.05. The condensed matrix therefore uses λ/s ≈ 0.05 μ/s². The direction predicts the limiting vertex will close by **4.0, 12.9, 10.3 and 12.4 times its gap**, and the slack bound cuts the whole box's step to **0.247, 0.077, 0.097 and 0.080**. The multiplier bound, the element bound and merit backtracking never bind. The coupled motion matches round 21: rows already pressed flat (g ≈ 3·10⁻⁴ d̂) are predicted to open by 16–19 times their gaps, while a far-side row closes.
   - **Where it does not match the hypothesis:** in frames 25–27 the weak stiffness is not the main cause. At the limiting vertex, IP's effective force μ/s is 1.7, 10.6 and 8.1 F_w, against a non-contact (momentum and elastic) push of 9.5, 17.0 and 17.8 F_w. At τ̂ = 1 the barrier level itself is too weak to stop the incoming momentum within the gap, whatever λ is.
   - **Frame 28 is the clean case:** the force is ample (49.6 against 8.7 F_w), and the stiffness deficit alone drives the closure.
2. **Did centering the multipliers change the limiting direction? Not at the selected frame (25); partly in frames 26–27; substantially only in frame 28.** At the same positions, rows, slacks and μ, with only λ_i := μ/s_i:

   | frame | carried: limiting motion / gap, slack bound | centered: same-vertex motion / gap, slack bound (limiting vertex) | IPC, same geometry: motion / gap, accepted step |
   |---|---|---|---|
   | **25 (selected)** | −4.02, **0.247** (v1) | −3.86, **0.258** (v1); the direction changes by 0.4 % | −0.14, 0.573 (CCD on v2) |
   | 26 | −12.9, **0.077** (v5) | −7.77, **0.128** (v5) | **+0.69**, 0.625 (CCD on v8) |
   | 27 | −10.3, **0.097** (v11) | −7.35, **0.135** (v11) | **+0.66**, 0.638 (CCD on v14) |
   | 28 | −12.4, **0.080** (v16) | −1.09, **0.512** (now v17) | **+0.80**, 1.0 |

   Even centered, the step stays slack-limited in all four frames. Replaying whole frames from identical states with central starts costs **more**, not less: 8, 10, 9 and 5 main systems against 7, 9, 8 and 6; over the 20 impact frames, 94 + 20 against 80 + 20.
3. **Did hard IP use more Newton solves than IPC? Yes, in impact and over the full run.**
   - **Impact (frames 24–43):** **80 main + 20 verification = 100 full-size solves**, against **67** for IPC.
   - **Full run (200 frames):** **415 + 200 = 615**, against **362**.
   - **How much the short-step mechanism costs:** in its frames (25–28), IP needs 30 main systems against IPC's 25 (+5). The other +8 impact main systems come from frames where every IP step is a full step: one system at τ̂ = 1 plus one or two at the final level, against IPC's 2–3. The verification systems (one per frame) are the largest single difference.

## Scene and shared numerics

| item | value |
|---|---|
| Body | 0.4 × 0.4 m elastic box, 16 × 16 cells, 289 vertices, 512 triangles (17 bottom vertices), alternating diagonals |
| Initial state | tilted 10° counter-clockwise, at rest, lowest vertex (the corner, v1) 0.3 m above the ground line y = 0; frictionless |
| Material | rest-stable 2D stable neo-Hookean Ψ = μ/2(‖F‖² − 2) − μ(J − 1) + λ/2(J − 1)²; E = 1e5 Pa, ν = 0.4 (plane-strain Lamé μ = 3.57e4, λ = 1.43e5); ρ = 1e3 kg/m³; lumped masses |
| Time | implicit Euler, dt = 0.01 s, 200 frames; x̃ = xₙ + dt vₙ + dt² g, g = −9.8 |
| Shared code | mesh, masses, `mech_eval` (energy, gradient, per-element PSD-projected Hessian), `elem_check` (exact first-inversion time; trials must keep every det > 0), `shared_solve` (sparse Cholesky), inertial prediction |
| Tolerances | d̂ = 0.01 m; tol_disp = velocity_tol·dt = 5e-4 m; force-balance target 5 F_w; ≤ 8 line-search trials; budget 256 main systems per frame |
| IPC κ | native scene-adaptive corridor (`global_contact_manager.cu`) evaluated in 2D: [1.1e6, 1.1e8]; the scene default 1e9 is clamped to κ = 1.1e8, κ·dt² = 1.1e4 |

**2D normalization.** The box is a plane-strain slab of unit depth, so masses, energies and forces are per metre of depth. The solver energy is ½‖x − x̃‖²_M + dt²Ψ_el (+ κdt²ΣB for IPC), and its gradient has the unit force·dt². **F_w** = mean vertex mass·|g|·dt² = 5.43e-4 is one vertex's weight in that unit; the whole body weighs 289 F_w.
- Force tolerances in F_w are per-vertex-weight units, so "5 F_w" means the same in 2D as in the 3D reference.
- The barrier has unit F_w·L with L = d̂: τ̂ = τ/(F_w L), τ₀ = 1, τ_final = tol_disp/L = 0.05.
- κ comes from the corridor and is proportional to the mean vertex mass, so κdt²/F_w does not depend on the depth convention.
- tol_disp is an absolute length for both solvers.

Everything resolved is in `results/settings.json`.

## Solvers

**IPC** (`src/ipc_solve_frame.m`):
- **Barrier:** the native log² barrier κdt²(D − d̂²)² log²(D/d̂²) with D = g² (`KappaBarrierLog2`), active for g < d̂. Its gradient and unprojected Hessian are evaluated at every iterate and every line-search trial.
- **Step:** half-plane CCD keeping (1 − 0.1) of each closing vertex's gap (`easy_vertex_half_plane_trajectory_filter.cu`), the shared element bound, and halving on the incremental energy (E ≤ E₀). CFL is off.
- **Stop (round 43 matched stop):** native |dx|∞ ≤ tol at the current state, plus at that state r_d ≤ 5 F_w with IPC's own law, full-correction CCD ≥ 1, valid elements and positive gaps. The state is returned without the terminal step.
- **Derivative check** (`results/derivative_checks_barrier.csv`): at gaps of 0.02 to 0.99 d̂, dB/dD, d²B/dD² and the assembled vertex gradient and Hessian agree with central differences to 1e-11–1e-8 (relative). The second difference of the energy agrees to 6e-7–1.6e-3; the largest error is at 0.99 d̂, where the stiffness is near zero. The normal stiffness is positive over the whole band.

**Hard IP, round-43 `base`** (`src/ip_solve_frame.m`; arm H + signed schedule + zero-CCD ground advance + early discovery):
- **Rows:** g_i − s_i = 0, s_iλ_i = μ, ∇E − Jᵀλ = 0, condensed to (H₀ + Jᵀdiag(λ/s)J)dx = −(∇E − Jᵀ(μ/s) + Jᵀdiag(λ/s)(g − s)).
- **Frame start:** τ̂ = 1, s = g, central λ = μ/s. Multipliers are carried from the preceding accepted frame by vertex identity (all 2,645 frame-start rows of the run were carried). The first barrier decrease goes directly to τ̂ = 0.05 when E_τ ≤ 10τ̂.
- **Step:** slack and multiplier fractions to the boundary (η = 0.995), update B, the split rule at τ̂ = 1 (5 uses), the element bound, s̃ = max(s + α ds, g(x̃)), and Armijo on F − τΣlog(s̃/d̂) + νΣ|g − s̃| with ν = 1.1 max|λ + dλ|.
- **Contacts:** no CCD. After every accepted step, all 289 vertex gaps are scanned (215 scans). In-band vertices become central rows (38 early rows). A nonpositive gap would trigger the halving ground advance; this never happened.
- **Stop:** the terminal update, then a fresh verification (one verification system): r_d with the actual λ, the correction, |g − s|, centrality, sign of s and λ, unresolved contacts, gaps and elements.
- **Not ported:** the element-redirect step (its trigger never armed; the element bound never limited a step), self contact, PCG (exact solves instead), and the rounding-certified determinant test (det > 0 instead). The reduced first step, P70 buffer and logarithmic gap update are left for future work, as asked.

## Matched accuracy (all 200 frames of both solvers qualify; no failed frame)

| at the returned state | hard IP | IPC |
|---|---|---|
| max ‖r_d‖∞ (own contact law; target 5 F_w) | 0.12 F_w (actual λ) | 4.91 F_w |
| max correction / tol_disp | 0.14 (verification system) | 0.91 (last main system) |
| lowest vertex gap | 2.8e-4 d̂ (row slacks equal gaps, max \|g − s\| = 0) | 0.51 d̂ |
| smallest det F / rest | 0.850 | 0.848 |
| IP only | max \|sλ/μ − 1\| = 1.19 (bound: 10τ_final, i.e. ≤ 10); unresolved contacts 0 | 20 native stops rejected by the matched check |

IP rests about 3 µm above the ground, while IPC rests at half the activation distance; this is the contact law, as in 3D. One IP system is counted as discarded: in frame 2 (free flight) the terminal correction was 1.6e-14 tol, the merit change was below rounding, all 8 trials were rejected, and the unchanged state passed verification.

## Newton solves (`results/phase_totals.csv`, `results/frames_*.csv`)

| phase | frames | IP main | IP verification | IP discarded / fallback | IP total | IPC main = total |
|---|---|---|---|---|---|---|
| free flight | 1–23 | 46 | 23 | 1 / 0 | 69 | 46 |
| impact | 24–43 | **80** | 20 | 0 / 0 | **100** | **67** |
| settling | 44–200 | 289 | 157 | 0 / 0 | 446 | 249 |
| **total** | 1–200 | **415** | 200 | 1 / 0 | **615** | **362** |

The phases follow the reference split: impact starts at the first frame in which either solver has a contact pair (frame 24 for both) and lasts 20 frames. Per-frame impact counts are in `results/summary.md`. The 3D reference showed the same ordering: 547 + 200 against 530 (matched stop), and in impact 107 + 20 against 94.

![solves per frame](results/plots/newton_solves_vs_frame.png)

## Same-state isolation (`results/replay_frame0NN*.{json,csv}`, `results/impact_sweep_first_limiting_iterates.csv`)

**Selection.** The rule was fixed before the run: the first impact frame with an accepted IP fraction ≤ 0.25. That is frame 25 (0.2456). Every impact frame was also replayed (a sweep) from IP's saved frame-start state (xₙ, vₙ, x̃) with IP's actual multiplier history. Frame 26 has the smallest fraction (0.077) and is written out in full as a supplement.
- **Replays:** IP reproduces its trajectory frames exactly (7 + 1 and 9 + 1 systems).
- **IPC from the same state:** 8 and 15 main systems, against 7 and 7 on its own trajectory. IP's flattened rows (g ≈ 5e-4 d̂) start IPC at r_d ≈ 2e6 F_w, so IPC replays after frame 25 are not like-for-like.

**Frame 25, first limiting iterate (k = 0, τ̂ = 1, one row).** The limiting vertex is v1, the corner, a carried row.

| quantity | value |
|---|---|
| g = s | 0.601 d̂ |
| λ | 0.0832 F_w |
| μ | 1 F_w L |
| sλ/μ | **0.050** |
| λ/s | 0.138 F_w/d̂ = **5 %** of μ/s² (2.77); 0.036 of the vertex mass |
| effective force μ/s | 1.66 F_w (non-contact push 9.48 F_w; IPC's force at the same gap 88 F_w, stiffness 704 F_w/d̂) |
| predicted normal motion | **−4.02 gaps** (closing) |

Restrictions recomputed from the raw direction:

| restriction | value |
|---|---|
| slack positivity | **0.2474** (0.995 × the raw crossing 0.2486, v1) |
| multiplier positivity | 1 (no row decreases) |
| element validity | none (no inversion on (0, 1]) |
| ground safety | exact crossing 0.2486, the same vertex; the early scan after the accept found every gap positive (no fallback) |
| merit backtracking | first trial accepted |

**Binding: slack positivity.**

**Assembly checks at s = g** (`assembly` in the JSON):
- K − H₀ is diagonal with λ_i/s_i on the rows' normal DOFs, to 1e-13.
- The right-hand-side contact term is −μ/s_i, to 2e-16.
- An unreduced (dx, ds, dλ) KKT solve reproduces the condensed direction to 2e-15.
- The physical force balance uses the actual λ: r_d = 37.6 F_w at frame 25. At frame 26, r_d = 39.3 F_w with λ, against 1,731 F_w if μ/s were used.

**Diagnostic work, counted apart from simulation:** 4 diagnostic solves per diagnosed frame (recomputed direction, KKT check, centered direction, IPC at the same geometry), 20 in the sweep. The same-state frame replays took 485 full-size solves.

![frame 25](results/plots/impact_frame_25_residual_and_step.png)

Supplementary frame 26: `results/plots/impact_frame_26_residual_and_step.png`.

## Files

- `run_experiment.m` is the single entry point.
- `src/` holds the code: `ipipc_params.m` (every setting), `scene_build.m`, `mech_eval.m`, `elem_check.m`, `shared_solve.m`, `ipc_barrier.m`, `ipc_contact_eval.m`, `ipc_solve_frame.m`, `ip_solve_frame.m`, `run_trajectory.m`, `check_derivatives.m`, `diagnose_impact.m`, `make_plots.m`.
- `results/` holds the outputs:
  - `settings.json`
  - `derivative_checks*.{csv,json}`
  - `frames_ip.csv`, `frames_ipc.csv` (per frame: counts, status, reason, accuracy)
  - `iters_ip.csv`, `iters_ipc.csv` (per iteration: τ̂, rows, r_d, step bounds, accepted fraction, limiter, limiting vertex, sλ/μ, closure, events)
  - `phase_totals.csv`
  - `impact_sweep_first_limiting_iterates.csv`
  - `replay_frame0{25,26,27,28,39}*.{json,csv}`
  - `summary.{md,json}`
  - `plots/*.png|.fig`
  - `trajectories.mat`
  - `run_log.txt`
