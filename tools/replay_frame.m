function [F, its, D] = replay_frame(source, frame, solver, diagnose)
%REPLAY_FRAME  Replay one saved frame from its frame-start state (positions, velocities, inertial prediction and, for
% IP, the multiplier history carried into that frame) and check it against the saved trajectory.
%   [F, its] = replay_frame('reference/baseline_2026-11-20/results', 26, 'ip')
%   [F, its] = replay_frame('runs/run_...', 26, 'ipc')
%   [F, its, D] = replay_frame(source, 26, 'ip', true)   % also: diagnose_impact.m (both replays, first limiting iterate)
% F is the frame record, its the per-iteration trace with captures (its(k).cap: x, rows, s, lambda, tau, dx, ds, ...).
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root, 'src'));
if nargin < 3, solver = 'ip'; end
if nargin < 4, diagnose = false; end
M = load(fullfile(char(source), 'trajectories.mat'), 'Rip', 'Ripc');
P = ipipc_params(); S = scene_build(P);
switch solver
    case 'ip'
        R = M.Rip; st = R.starts{frame};
        [~, F, its] = ip_solve_frame(S, P, st.x, st.xt, frame, st.hist, true);
    case 'ipc'
        R = M.Ripc; st = R.starts{frame};
        [~, F, its] = ipc_solve_frame(S, P, st.x, st.xt, frame, true);
end
ref = R.iters([R.iters.frame] == frame);
same = numel(ref) == numel(its) && isequaln([ref.alpha_eff], [its.alpha_eff]) && isequaln([ref.rd_inf], [its.rd_inf]);
fprintf('%s frame %d: replay %d main + %d verification (%s); saved trajectory %d + %d (%s); iterations identical: %d\n', ...
    solver, frame, F.main, F.verification, F.status, R.frames(frame).main, R.frames(frame).verification, ...
    R.frames(frame).status, same);
D = [];
if diagnose
    D = diagnose_impact(S, P, M.Rip, M.Ripc, frame, 'replay_frame');
    if D.limiting.found
        fprintf('first limiting IP iterate k=%d: %s, accepted %.4g\n', D.limiting.k, D.limiting.limiter, D.limiting.carried.accepted_alpha_eff);
    else
        fprintf('%s\n', D.limiting.reason);
    end
end
end
