# Auto-generated summary (run_experiment.m)

Phases: free flight = frames before the first frame in which either solver has a contact pair; impact = that frame and the next 19; settling = the rest. First contact frame 24.

Selected frame: 25 (first impact frame with an accepted IP position fraction <= 0.25 (it is 0.2456)).

## Newton solves by phase

| phase | frames | solver | main | verification | discarded | fallback | total full-size | failed frames |
|---|---|---|---|---|---|---|---|---|
| free flight | 1-23 | ip | 46 | 23 | 1 | 0 | 69 | 0  |
| free flight | 1-23 | ipc | 46 | 0 | 0 | 0 | 46 | 0  |
| impact | 24-43 | ip | 80 | 20 | 0 | 0 | 100 | 0  |
| impact | 24-43 | ipc | 67 | 0 | 0 | 0 | 67 | 0  |
| settling | 44-200 | ip | 289 | 157 | 0 | 0 | 446 | 0  |
| settling | 44-200 | ipc | 249 | 0 | 0 | 0 | 249 | 0  |
| total | 1-200 | ip | 415 | 200 | 1 | 0 | 615 | 0  |
| total | 1-200 | ipc | 362 | 0 | 0 | 0 | 362 | 0  |

## Impact frames, per frame

| frame | IP main | IP ver | IP min accepted fraction | IPC main | IPC min accepted fraction |
|---|---|---|---|---|---|
| 24 | 3 | 1 | 1 | 4 | 1 |
| 25 | 7 | 1 | 0.2456 | 7 | 0.5659 |
| 26 | 9 | 1 | 0.07715 | 7 | 0.5876 |
| 27 | 8 | 1 | 0.09656 | 8 | 0.5551 |
| 28 | 6 | 1 | 0.08003 | 3 | 1 |
| 29 | 3 | 1 | 1 | 3 | 1 |
| 30 | 3 | 1 | 1 | 3 | 1 |
| 31 | 3 | 1 | 1 | 2 | 1 |
| 32 | 3 | 1 | 1 | 2 | 1 |
| 33 | 3 | 1 | 1 | 3 | 1 |
| 34 | 4 | 1 | 1 | 4 | 1 |
| 35 | 4 | 1 | 1 | 3 | 1 |
| 36 | 3 | 1 | 1 | 3 | 1 |
| 37 | 4 | 1 | 1 | 3 | 1 |
| 38 | 4 | 1 | 1 | 2 | 1 |
| 39 | 4 | 1 | 0.9808 | 2 | 1 |
| 40 | 3 | 1 | 1 | 2 | 1 |
| 41 | 2 | 1 | 1 | 2 | 1 |
| 42 | 2 | 1 | 1 | 2 | 1 |
| 43 | 2 | 1 | 1 | 2 | 1 |

## Accuracy at the returned states

- max ||r_d||_inf: IP 0.1218 F_w, IPC 4.906 F_w (target 5)
- max returned correction / tol: IP 0.139 (verification system), IPC 0.9085 (last main system)
- min vertex gap: IP 0.0002827 d_hat, IPC 0.5135 d_hat; min det F (rel. rest): IP 0.8504, IPC 0.8476
- IP: max |s lambda/tau - 1| 1.189, max |g - s| 0 d_hat, max unresolved contacts 0
- failed frames: IP [], IPC []

## Same-state replay of frame 25 (selected)

- IP replay 7 + 1 (trajectory 7 + 1, same iterations 1, max diff alpha 0); IPC replay 8 main (own trajectory 7); IP central start 8 + 1
- first limiting iterate k=0, tau_hat 1, 1 rows; accepted 0.2474 (slack positivity)
- vertex 1 (carried (history)): g = 0.6013 d_hat, s = 0.6013 d_hat, lambda = 0.08315 F_w, mu = 1 F_w L, s lambda/mu = 0.05, lambda/s = 0.1383 F_w/d_hat (central mu/s^2 = 2.765), predicted motion -4.022 x gap
- restrictions: slack 0.2474 (v1), multiplier 1 (vNaN), element Inf, ground exact crossing 0.2486 (v1), merit accepted 0.2474 after 0 backtracks; binding: slack positivity
- centered lambda = mu/s: slack 0.2578 (v1), same-vertex motion -3.859 x gap, max closure 3.859 x gap, merit accepted 0.2578, binding slack positivity
- IPC same geometry (1 active, same pairs 1): same-vertex motion -0.1389 x gap, CCD 0.573 (v2), accepted 0.573
- assembly: |g - s| 0 d_hat, stiffness rel err 1.07e-13, rhs mu/s rel err 2.4e-16, off-diag 0, KKT dx/ds/dlambda rel err 1.74e-15 / 2.73e-15 / 6.26e-16; r_d with actual lambda 37.63 F_w vs with mu/s 37.63 F_w
- diagnostic solves 4, replay solves 25 (not simulation work)

## Same-state replay of frame 26 (impact frame with the smallest accepted IP fraction)

- IP replay 9 + 1 (trajectory 9 + 1, same iterations 1, max diff alpha 0); IPC replay 15 main (own trajectory 7); IP central start 10 + 1
- first limiting iterate k=0, tau_hat 1, 7 rows; accepted 0.07715 (slack positivity)
- vertex 5 (carried (history)): g = 0.09471 d_hat, s = 0.09471 d_hat, lambda = 0.5923 F_w, mu = 1 F_w L, s lambda/mu = 0.0561, lambda/s = 6.254 F_w/d_hat (central mu/s^2 = 111.5), predicted motion -12.9 x gap
- restrictions: slack 0.07715 (v5), multiplier 1 (vNaN), element Inf, ground exact crossing 0.07754 (v5), merit accepted 0.07715 after 0 backtracks; binding: slack positivity
- centered lambda = mu/s: slack 0.128 (v5), same-vertex motion -7.774 x gap, max closure 7.774 x gap, merit accepted 0.128, binding slack positivity
- IPC same geometry (7 active, same pairs 1): same-vertex motion 0.6923 x gap, CCD 0.6245 (v8), accepted 0.6245
- assembly: |g - s| 0 d_hat, stiffness rel err 8.21e-14, rhs mu/s rel err 1.71e-16, off-diag 0, KKT dx/ds/dlambda rel err 2.07e-15 / 1.58e-15 / 3.33e-15; r_d with actual lambda 39.31 F_w vs with mu/s 1731 F_w
- diagnostic solves 4, replay solves 36 (not simulation work)

## Impact sweep: every impact frame replayed from IP's frame-start state

| frame | IP main+ver | IP central start | IPC from IP state | IPC own | k | tau_hat | vertex | source | s lam/mu | motion/gap carried | slack bound carried | binding | motion/gap centered | slack bound centered | motion/gap IPC |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 24 | 3+1 | 3+1 | 4 | 4 | NaN | NaN | NaN |  | NaN | NaN | NaN | none: full steps | NaN | NaN | NaN |
| 25 | 7+1 | 8+1 | 8 | 7 | 0 | 1 | 1 | carried (history) | 0.05 | -4.02 | 0.247 | slack positivity | -3.86 | 0.258 | -0.139 |
| 26 | 9+1 | 10+1 | 15 | 7 | 0 | 1 | 5 | carried (history) | 0.0561 | -12.9 | 0.0771 | slack positivity | -7.77 | 0.128 | 0.692 |
| 27 | 8+1 | 9+1 | 16 | 8 | 0 | 1 | 11 | carried (history) | 0.0504 | -10.3 | 0.0966 | slack positivity | -7.35 | 0.135 | 0.656 |
| 28 | 6+1 | 5+1 | 16 | 3 | 0 | 1 | 16 | carried (history) | 0.0499 | -12.4 | 0.08 | slack positivity | -1.09 | 0.512 | 0.797 |
| 29 | 3+1 | 3+1 | 16 | 3 | NaN | NaN | NaN |  | NaN | NaN | NaN | none: full steps | NaN | NaN | NaN |
| 30 | 3+1 | 3+1 | 16 | 3 | NaN | NaN | NaN |  | NaN | NaN | NaN | none: full steps | NaN | NaN | NaN |
| 31 | 3+1 | 2+1 | 16 | 2 | NaN | NaN | NaN |  | NaN | NaN | NaN | none: full steps | NaN | NaN | NaN |
| 32 | 3+1 | 2+1 | 15 | 2 | NaN | NaN | NaN |  | NaN | NaN | NaN | none: full steps | NaN | NaN | NaN |
| 33 | 3+1 | 2+1 | 15 | 3 | NaN | NaN | NaN |  | NaN | NaN | NaN | none: full steps | NaN | NaN | NaN |
| 34 | 4+1 | 6+1 | 14 | 4 | NaN | NaN | NaN |  | NaN | NaN | NaN | none: full steps | NaN | NaN | NaN |
| 35 | 4+1 | 7+1 | 14 | 3 | NaN | NaN | NaN |  | NaN | NaN | NaN | none: full steps | NaN | NaN | NaN |
| 36 | 3+1 | 6+1 | 14 | 3 | NaN | NaN | NaN |  | NaN | NaN | NaN | none: full steps | NaN | NaN | NaN |
| 37 | 4+1 | 6+1 | 13 | 3 | NaN | NaN | NaN |  | NaN | NaN | NaN | none: full steps | NaN | NaN | NaN |
| 38 | 4+1 | 6+1 | 13 | 2 | NaN | NaN | NaN |  | NaN | NaN | NaN | none: full steps | NaN | NaN | NaN |
| 39 | 4+1 | 4+1 | 13 | 2 | 1 | 0.05 | 17 | carried (history) | 16.9 | -1.01 | 0.981 | slack positivity | -5.23 | 0.19 | 0.756 |
| 40 | 3+1 | 3+1 | 13 | 2 | NaN | NaN | NaN |  | NaN | NaN | NaN | none: full steps | NaN | NaN | NaN |
| 41 | 2+1 | 3+1 | 13 | 2 | NaN | NaN | NaN |  | NaN | NaN | NaN | none: full steps | NaN | NaN | NaN |
| 42 | 2+1 | 3+1 | 13 | 2 | NaN | NaN | NaN |  | NaN | NaN | NaN | none: full steps | NaN | NaN | NaN |
| 43 | 2+1 | 3+1 | 14 | 2 | NaN | NaN | NaN |  | NaN | NaN | NaN | none: full steps | NaN | NaN | NaN |

Totals over the 20 impact frames: IP 80 + 20, IP central start 94 + 20, IPC from IP's states 271, IPC on its own trajectory 67. Diagnostic solves 20, replay solves 485.

## Answers (rule-based)

1. Short-step mechanism reproduced: 1. frame 25, iterate k=0 (tau_hat 1): binding restriction slack positivity at vertex 1 (carried (history)); s*lambda/mu = 0.05, lambda/s = 0.05 of the central mu/s^2; predicted normal motion -4.02 x its gap; accepted fraction 0.2474
2. Centering changed the limiting direction: 0. carried: slack bound 0.2474 (vertex 1), closure 4.02 x gap; centered (lambda = mu/s): slack bound 0.2578 (vertex 1), same vertex motion -3.86 x gap, max closure 3.86 x gap, direction change 0.00391 (rel.); IPC same geometry: motion -0.139 x gap, accepted 0.573
3. Impact: IP 80 main + 20 ver = 100 vs IPC 67. Full run: IP 415 main + 200 ver = 615 vs IPC 362.

Plots: E:\ZiyanXiong\github\ipipc-matlab\results\plots\newton_solves_vs_frame.png, E:\ZiyanXiong\github\ipipc-matlab\results\plots\impact_frame_25_residual_and_step.png, E:\ZiyanXiong\github\ipipc-matlab\results\plots\impact_frame_26_residual_and_step.png
