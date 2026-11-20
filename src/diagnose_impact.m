function D = diagnose_impact(S, P, Rip, Ripc, f, label)
%DIAGNOSE_IMPACT  Same-state replays of frame f and the analysis of its first limiting IP iterate.
% Every replay starts from a saved frame-start state (positions, velocities, inertial prediction):
%   ip_replay          IP at IP's own frame-start state, with its actual multiplier history -> agreement with its trajectory
%   ipc_replay         IPC at IPC's own frame-start state -> agreement with IPC's trajectory
%   ipc_at_ip_state    IPC at IP's frame-start state (a cross replay at the same state, not a trajectory)
%   ip_central_restart IP at IP's frame-start state with every multiplier started centrally (history ignored)
% The first limiting iterate is the first main IP iterate of the replay whose accepted position fraction is in (0, 1);
% diagnose_iterate.m analyses it whatever restriction set it. Replay and diagnostic solves are counted, not simulation work.
if nargin < 6, label = 'impact sweep'; end
assert(f <= numel(Rip.starts), 'diagnose_impact: IP did not reach frame %d', f);
st = Rip.starts{f};
D = struct('label', label, 'frame', f, 'diagnostic_solves', 0, 'replay_solves', 0);

[~, Fi, iti] = ip_solve_frame(S, P, st.x, st.xt, f, st.hist, true);
D.ip_replay = agreement('ip', f, Fi, iti, Rip);
if f <= numel(Ripc.starts)
    stc = Ripc.starts{f};
    [~, Fs, its_] = ipc_solve_frame(S, P, stc.x, stc.xt, f, false);
    D.ipc_replay = agreement('ipc', f, Fs, its_, Ripc);
    n_self = Fs.main;
else
    D.ipc_replay = agreement('ipc', f, [], [], Ripc); n_self = 0;
end
[~, Fc, itc] = ipc_solve_frame(S, P, st.x, st.xt, f, true);
own = NaN; if f <= numel(Ripc.frames), own = Ripc.frames(f).main; end
D.ipc_at_ip_state = struct('main', Fc.main, 'status', Fc.status, 'ipc_own_trajectory_main', own);
[~, Fz] = ip_solve_frame(S, P, st.x, st.xt, f, [], false);
D.ip_central_restart = struct('main', Fz.main, 'verification', Fz.verification, 'status', Fz.status);
D.same_state = struct('max_abs_diff_first_iterate_x', max(abs(iti(1).cap.x - itc(1).cap.x)), ...
    'max_abs_diff_xt_vs_prediction', max(abs(iti(1).cap.xt - (st.x + S.dt * st.v + S.dt2 * S.gvec))), ...
    'hist_rows', numel(field_or(st.hist, 'rv', [])), 'hist_frame', field_or(st.hist, 'frame', NaN), ...
    'rows_carried_at_start', Fi.hist_reused, 'rows_central_at_start', Fi.hist_fresh);
D.replay_solves = Fi.main + Fi.verification + n_self + Fc.main + Fz.main + Fz.verification;
D.ip_iters = strip_cap(iti); D.ipc_iters = strip_cap(itc);

mains = iti(strcmp({iti.kind}, 'main'));
q = find([mains.alpha_eff] < 1 & [mains.alpha_eff] > 0, 1);
if isempty(q)
    D.limiting = struct('found', false, 'reason', 'every accepted IP step of this frame is the full step');
    return;
end
D.limiting = diagnose_iterate(S, P, mains(q));
D.diagnostic_solves = D.limiting.diagnostic_solves;
end

% ======================================================================================== helpers
function A = agreement(solver, f, Fr, itr, R)
% the replay against the same solver's own trajectory frame (exact equality expected: deterministic code)
A = struct('solver', solver, 'frame', f, 'available', false, 'main_replay', NaN, 'main_trajectory', NaN, ...
    'verification_replay', NaN, 'verification_trajectory', NaN, 'status_replay', '', 'status_trajectory', '', ...
    'same_iteration_count', false, 'max_abs_diff_alpha', NaN, 'max_abs_diff_alpha_eff', NaN, ...
    'max_abs_diff_rd_inf', NaN, 'max_abs_diff_dx_tol', NaN, 'max_abs_diff_min_gap', NaN, 'same_limiters', false, ...
    'identical', false);
if isempty(Fr) || f > numel(R.frames), return; end
ref = R.iters([R.iters.frame] == f);
A.available = true;
A.main_replay = Fr.main; A.main_trajectory = R.frames(f).main;
A.verification_replay = Fr.verification; A.verification_trajectory = R.frames(f).verification;
A.status_replay = Fr.status; A.status_trajectory = R.frames(f).status;
A.same_iteration_count = numel(itr) == numel(ref);
if A.same_iteration_count
    d = @(name) max([0, abs(z([itr.(name)]) - z([ref.(name)]))]);
    A.max_abs_diff_alpha = d('alpha'); A.max_abs_diff_alpha_eff = d('alpha_eff');
    A.max_abs_diff_rd_inf = d('rd_inf'); A.max_abs_diff_dx_tol = d('dx_tol'); A.max_abs_diff_min_gap = d('min_gap_dhat');
    A.same_limiters = isequal({itr.limiter}, {ref.limiter});
end
A.identical = A.same_iteration_count && A.main_replay == A.main_trajectory && ...
    A.verification_replay == A.verification_trajectory && strcmp(A.status_replay, A.status_trajectory) && ...
    A.max_abs_diff_alpha == 0 && A.max_abs_diff_alpha_eff == 0 && A.max_abs_diff_rd_inf == 0 && ...
    A.max_abs_diff_dx_tol == 0 && A.max_abs_diff_min_gap == 0 && A.same_limiters;
end

function v = z(v)
v(isnan(v)) = 0; v(isinf(v)) = realmax;
end

function its = strip_cap(its)
for q = 1:numel(its), its(q).cap = []; end
end

function v = field_or(s, f, d)
if isempty(s) || ~isfield(s, f), v = d; else, v = s.(f); end
end
