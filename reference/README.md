# Preserved reference evidence

Nothing in this folder is written by the code. `run_experiment` refuses any output path inside `reference/`.

## `baseline_2026-11-20/` — the verified 2D run

| path | contents |
|---|---|
| `results/` | the complete verified run, byte-identical to the original `results/` folder: `trajectories.mat` (full MATLAB state, including every frame-start state used for replays), resolved `settings.json`, per-frame and per-iteration traces, the impact-frame replay diagnoses, plots, `summary.{md,json}`, `run_log.txt` |
| `code/` | the exact source that produced `results/` (`run_experiment.m`, `src/`), and the report as first written (`REPORT.md`, since superseded by `../../REPORT.md`) |

Some text files inside `results/` (`summary.md`, `run_log.txt`) mention the folder's original absolute location `E:\ZiyanXiong\github\ipipc-matlab\results\`. The folder was moved here without any content change.

**Integrity.** `baseline_2026-11-20.sha256` lists all 54 files (39 in `results/`, 15 in `code/`) in `sha256sum` format. The manifest file itself has SHA-256 `28759a53183f73ff3966ef7ac1e2dc2b31d201725d9d684346a5061f83d1b96e`. To verify:

```
cd reference && sha256sum -c baseline_2026-11-20.sha256      # Git Bash / Linux / macOS
```

```matlab
addpath tools; verify_reference                             % MATLAB, from the project folder
```

**Expected values.** 200/200 qualified frames for both solvers. Hard IP: 415 main + 200 verification = 615 full-size solves. IPC: 362 main = 362. Impact frames 24–43: IP 80 + 20, IPC 67. `tools/compare_to_reference.m` checks a new run against everything in `results/`.

## `libuipc_ip_sources.*` — the 3D reference the port was read from

The 2D port was read from `E:\ZiyanXiong\github\libuipc-ip` (read-only) on 2026-11-20. That checkout is at commit `b611ae72ead0228b8d49cec7703c42bd43a1d40f` on branch `lagged-friction`, but **the commit alone does not identify the sources used**:
- **untracked:** `ip_solver.cu`, `ip_logic.h` and `ip_multiplier_history.h` (the hard-IP solver) and ADR 0013;
- **uncommitted edits:** `advance_ipc.cu`, `ipc_vertex_half_plane_contact_function.h`, `scene_default_config.cpp`, `eval_scenes.py` and `presets.py`;
- **git-ignored:** the round reports and runners under `local/`.

`libuipc_ip_sources.sha256` records the SHA-256 of every file the port was read from, with paths relative to the libuipc-ip root. `src/backends/cuda/engine/ip_solver.cu` hashes to `a58e02d6…412e1`, the same digest that `local/ip_vs_ipc_block_drop/round43/identity.txt` records for the round-43 build. The hard-IP logic was therefore read from the solver source that produced round 43's `base` arm.

| role | file (libuipc-ip) | git status on 2026-11-20 |
|---|---|---|
| hard IP, one frame | `src/backends/cuda/engine/ip_solver.cu` | untracked |
| barrier schedule, merit, Armijo, fractions to the boundary | `src/backends/cuda/engine/ip_logic.h` | untracked |
| multiplier history | `src/backends/cuda/engine/ip_multiplier_history.h` | untracked |
| native IPC Newton loop and matched stop | `src/backends/cuda/engine/advance_ipc.cu` | modified |
| κ corridor | `src/backends/cuda/contact_system/global_contact_manager.cu` | clean |
| log² barrier and derivatives | `src/backends/cuda/contact_system/contact_models/sym/codim_ipc_contact.inl` | clean |
| vertex–half-plane contact | `.../contact_models/ipc_vertex_half_plane_contact_function.h` | modified |
| half-plane distance | `.../contact_models/sym/vertex_half_plane_distance.inl` | clean |
| half-plane CCD (η = 0.1) | `src/backends/cuda/collision_detection/filters/easy_vertex_half_plane_trajectory_filter.cu` | clean |
| scene diagonal | `src/backends/cuda/global_geometry/global_vertex_manager.cu` | clean |
| defaults (d̂, tolerances, line search, CCD tol, CFL) | `src/core/core/scene_default_config.cpp` | modified |
| 3D scene `fem_block_drop` | `benchmarks/converged/eval_scenes.py` | modified |
| presets `ip-matched`, `ipc-matched` | `benchmarks/converged/presets.py` | modified |
| maintained-solver decision | `agent_docs/adr/0013-maintained-ip-solver.md` | untracked |
| runs and arms (H, `base`) | `local/ip_vs_ipc_block_drop/run.sh`, `round43/run43s.sh` | ignored (`local/`) |
| evidence | `local/ip_vs_ipc_block_drop/round20/REPORT.md`, `round21/REPORT.md`, `round43/REPORT.md`, `round43/identity.txt` | ignored (`local/`) |

The round-43 `base` arm is `run.sh` arm `H` (preset `ip-matched`, `UIPC_IP_MULTIPLIER_HISTORY=1`, `UIPC_IP_TERMINAL_UPDATE=1`) plus `UIPC_IP_SIGNED_OUTER=1`, `UIPC_IP_SIGNED_OUTER_ADVANCE=ground` and `UIPC_IP_SIGNED_EARLY_DISCOVERY=1`.
