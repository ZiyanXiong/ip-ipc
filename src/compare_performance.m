function C = compare_performance(Rip, Ripc, ph)
%COMPARE_PERFORMANCE  Solve-count and wall-time comparison of two trajectories, only when BOTH are complete.
% A trajectory that stopped at a failed frame has fewer frames, so its totals are smaller by construction. Such a run
% is never compared: C.comparable = false, every comparison field is NaN/empty, and C.verdict says why. The per-phase
% table still lists each solver's own counts over its qualified frames, flagged as incomplete.
C = struct();
C.ip = status_of(Rip); C.ipc = status_of(Ripc);
C.comparable = C.ip.complete && C.ipc.complete;
names = {'free flight', 'impact', 'settling', 'total'};
sets = {ph.free, ph.impact, ph.settling, ph.all};
rows = {};
for q = 1:numel(names)
    fs = sets{q};
    if isempty(fs), continue; end
    for R = {Rip, Ripc}
        r = R{1};
        have = intersect(fs, [r.frames([r.frames.qualified]).frame]);
        F = r.frames(ismember([r.frames.frame], have));
        failed = [r.frames(~[r.frames.qualified] & ismember([r.frames.frame], fs)).frame];
        missing = setdiff(fs, [r.frames.frame]);
        rows(end+1, :) = {names{q}, sprintf('%d-%d', fs(1), fs(end)), r.solver, numel(fs), numel(have), ...
            r.complete, sum([F.main]), sum([F.verification]), sum([F.discarded]), sum([F.fallback]), ...
            sum([F.main]) + sum([F.verification]), numel(failed), list(failed), numel(missing), sum([F.wall_s])}; %#ok<AGROW>
    end
end
C.table = cell2table(rows, 'VariableNames', {'phase', 'frames', 'solver', 'n_frames', 'qualified_frames', ...
    'trajectory_complete', 'main_systems', 'verification_systems', 'discarded_systems', 'fallback_systems', ...
    'total_full_size_solves', 'failed_frames', 'failed_frame_list', 'frames_not_run', 'wall_s'});
C.phases = struct();
for q = 1:numel(names)
    k = matlab.lang.makeValidName(names{q});
    t = C.table(strcmp(C.table.phase, names{q}), :);
    if isempty(t) || ~C.comparable
        C.phases.(k) = struct('ip_main', NaN, 'ip_verification', NaN, 'ip_total', NaN, 'ipc_main', NaN, 'ipc_total', NaN, ...
            'difference_main', NaN, 'difference_total', NaN, 'ip_wall_s', NaN, 'ipc_wall_s', NaN);
        continue;
    end
    a = t(strcmp(t.solver, 'ip'), :); b = t(strcmp(t.solver, 'ipc'), :);
    C.phases.(k) = struct('ip_main', a.main_systems, 'ip_verification', a.verification_systems, ...
        'ip_total', a.total_full_size_solves, 'ipc_main', b.main_systems, 'ipc_total', b.total_full_size_solves, ...
        'difference_main', a.main_systems - b.main_systems, 'difference_total', a.total_full_size_solves - b.total_full_size_solves, ...
        'ip_wall_s', a.wall_s, 'ipc_wall_s', b.wall_s);
end
if C.comparable
    i = C.phases.impact; t = C.phases.total;
    C.verdict = sprintf(['both trajectories complete (%d/%d frames qualified each). Impact: IP %d main + %d verification = %d, ' ...
        'IPC %d. Full run: IP %d + %d = %d, IPC %d.'], C.ip.last_qualified_frame, C.ip.frames_requested, i.ip_main, ...
        i.ip_verification, i.ip_total, i.ipc_main, t.ip_main, t.ip_verification, t.ip_total, t.ipc_main);
else
    why = {};
    if ~C.ip.complete, why{end+1} = sprintf('IP incomplete: first failure at frame %d (%s)', C.ip.first_failure, C.ip.failure_reason); end
    if ~C.ipc.complete, why{end+1} = sprintf('IPC incomplete: first failure at frame %d (%s)', C.ipc.first_failure, C.ipc.failure_reason); end
    C.verdict = ['not comparable: ' strjoin(why, '; ') '. No solve-count or wall-time comparison is reported for a partial run.'];
end
end

function s = status_of(R)
s = struct('complete', R.complete, 'frames_requested', R.frames_requested, 'last_qualified_frame', R.last_qualified_frame, ...
    'first_failure', R.first_failure, 'failure_reason', '');
if ~isempty(R.failure), s.failure_reason = sprintf('%s: %s', R.failure.status, R.failure.reason); end
end

function t = list(v)
t = strjoin(arrayfun(@(z) sprintf('%d', z), v, 'UniformOutput', false), ' ');
end
