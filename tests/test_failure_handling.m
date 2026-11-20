function T = test_failure_handling()
%TEST_FAILURE_HANDLING  A controlled failure exercises the failed-frame path (run_trajectory.m, compare_performance.m).
% The fault gives one frame a budget of 1 main system (the solver's own budget stop), nothing else changes:
%   IP  fails at frame 25 (impact; needs 7 + 1 systems), IPC fails at frame 5 (free flight; needs 2).
% Checks: the failing frame is identified; its diagnostics are kept; the last qualified state (positions, velocity,
% multiplier history, time) is preserved and equals an unfaulted run stopped there; the state is not advanced; and a
% partial run, which has fewer solves and less wall time by construction, is never reported as a faster completion.
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root, 'src'));
P = ipipc_params(); S = scene_build(P);
nf = 40; P40 = P; P40.frames = nf;
T = table('Size', [0 3], 'VariableTypes', {'string', 'logical', 'string'}, 'VariableNames', {'check', 'passed', 'detail'});

% ---- IP, controlled failure at frame 25
Rf = run_trajectory(S, P, 'ip', nf, struct('fault', struct('frame', 25, 'cap_main', 1)));
q = [Rf.frames.qualified];
T = add(T, "IP: failing frame identified", ~Rf.complete && Rf.first_failure == 25 && Rf.last_qualified_frame == 24 && ...
    numel(Rf.frames) == 25 && all(q(1:24)) && ~q(25), sprintf('first failure %d, last qualified %d, frames attempted %d', ...
    Rf.first_failure, Rf.last_qualified_frame, numel(Rf.frames)));
T = add(T, "IP: failure diagnostics preserved", ~isempty(Rf.failure) && strcmp(Rf.failure.status, 'budget') && ...
    Rf.failure.fault_injected && Rf.frames(25).main == 1 && ~isempty(Rf.failure.iters) && ...
    all([Rf.failure.iters.frame] == 25) && numel(Rf.failure.returned_state) == S.ndof, ...
    sprintf('status %s (%s); %d iterations kept', Rf.failure.status, Rf.failure.reason, numel(Rf.failure.iters)));
st = Rf.starts{25};
T = add(T, "IP: last qualified state preserved, not advanced", isequal(Rf.x_end, st.x) && isequal(Rf.v_end, st.v) && ...
    isequal(Rf.hist_end, st.hist) && abs(Rf.time_end - 24 * S.dt) < 1e-15 && ~Rf.frames(25).state_advanced && ...
    ~isequal(Rf.failure.returned_state, Rf.x_end), sprintf('time_end %.2f s, history rows %d, failed attempt moved %.3g m', ...
    Rf.time_end, numel(Rf.hist_end.rv), max(abs(Rf.failure.returned_state - Rf.x_end))));
Rc = run_trajectory(S, P, 'ip', 24);
T = add(T, "IP: equals an unfaulted run stopped at frame 24", Rc.complete && isequal(Rc.x_end, Rf.x_end) && ...
    isequal(Rc.v_end, Rf.v_end) && isequal(Rc.hist_end, Rf.hist_end), 'x, v and multiplier history bitwise equal');

% ---- the comparison refuses the partial run
Ripc = run_trajectory(S, P, 'ipc', nf);
ph = trajectory_phases(Rf, Ripc, P40);
C = compare_performance(Rf, Ripc, ph);
naive_ip = sum([Rf.frames(q).total]); naive_ipc = sum([Ripc.frames.total]);
wall_ip = sum([Rf.frames.wall_s]); wall_ipc = sum([Ripc.frames.wall_s]);
T = add(T, "partial IP run has fewer solves than the complete IPC run", naive_ip < naive_ipc, ...
    sprintf(['IP %d solves over its 24 qualified frames vs IPC %d over %d frames (a naive total would call IP cheaper); ' ...
    'wall time %.2f s vs %.2f s (information only)'], naive_ip, naive_ipc, nf, wall_ip, wall_ipc));
T = add(T, "partial run is not reported as a faster completion", ~C.comparable && contains(C.verdict, 'not comparable') && ...
    ~contains(lower(C.verdict), 'faster') && isnan(C.phases.total.ip_total) && isnan(C.phases.total.difference_total) && ...
    isnan(C.phases.total.ip_wall_s) && any(~C.table.trajectory_complete) && ...
    all(C.table.frames_not_run(strcmp(C.table.solver, 'ip') & strcmp(C.table.phase, 'total')) == nf - 25), C.verdict);

% ---- IPC, controlled failure at frame 5
Rg = run_trajectory(S, P, 'ipc', nf, struct('fault', struct('frame', 5, 'cap_main', 1)));
st = Rg.starts{5};
T = add(T, "IPC: failing frame identified, state preserved", ~Rg.complete && Rg.first_failure == 5 && ...
    Rg.last_qualified_frame == 4 && isequal(Rg.x_end, st.x) && isequal(Rg.v_end, st.v) && strcmp(Rg.failure.status, 'budget'), ...
    sprintf('first failure %d (%s)', Rg.first_failure, Rg.failure.reason));
Rip = run_trajectory(S, P, 'ip', nf);
C2 = compare_performance(Rip, Rg, trajectory_phases(Rip, Rg, P40));
T = add(T, "partial IPC run is not compared either", ~C2.comparable && contains(C2.verdict, 'IPC incomplete'), C2.verdict);

% ---- both complete: comparable
C3 = compare_performance(Rip, Ripc, trajectory_phases(Rip, Ripc, P40));
T = add(T, "two complete runs are compared", C3.comparable && isfinite(C3.phases.total.ip_total), C3.verdict);

disp(T);
if ~all(T.passed), error('test_failure_handling: %d of %d checks failed', sum(~T.passed), height(T)); end
fprintf('test_failure_handling: all %d checks passed\n', height(T));
end

function T = add(T, name, ok, detail)
T = [T; {name, logical(ok), string(detail)}];
end
