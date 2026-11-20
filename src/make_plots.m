function files = make_plots(Rip, Ripc, ph, Dsel, outdir, Dextra)
%MAKE_PLOTS  (1) Newton solves per frame for both solvers; (2) the selected impact frame (same-state replay):
% force residual and accepted step fraction versus solve number. Static PNG + FIG. Colors: reference palette
% categorical slots 1-2 (IP blue #2a78d6, IPC orange #eb6834); text in neutral ink; phases as gray bands.
% Each solver is drawn over the frames it attempted; a failed (unqualified) frame is marked with x and ends its line.
% Totals in the legend are over qualified frames; "(incomplete)" marks a trajectory that stopped at a failure.
if ~exist(outdir, 'dir'), mkdir(outdir); end
cIP = hex('#2a78d6'); cIPC = hex('#eb6834');
ink = hex('#0b0b0b'); ink2 = hex('#52514e'); band = hex('#efeeea');
files = {};

% ------------------------------------------------------------------ 1. solves per frame
fi = [Rip.frames.frame]; fc = [Ripc.frames.frame];
ip_tot = [Rip.frames.total]; ip_main = [Rip.frames.main]; ipc_main = [Ripc.frames.main];
qi = [Rip.frames.qualified]; qc = [Ripc.frames.qualified];
nf = max([ph.all, fi, fc]);
tag = @(R) ternary(R.complete, '', ' (incomplete)');
f1 = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1100 720]);
tl1 = tiledlayout(f1, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
ax = nexttile(tl1); hold(ax, 'on');
ymax = max([ip_tot, ipc_main, 1]) + 2;
if ~isempty(ph.impact), draw_band(ax, ph.impact, ymax, band); end
stairs(ax, [fi fi(end)+1] - 0.5, [ip_main ip_main(end)], '--', 'Color', cIP, 'LineWidth', 1);
stairs(ax, [fc fc(end)+1] - 0.5, [ipc_main ipc_main(end)], '-', 'Color', cIPC, 'LineWidth', 2);
stairs(ax, [fi fi(end)+1] - 0.5, [ip_tot ip_tot(end)], '-', 'Color', cIP, 'LineWidth', 2);
if any(~qi), plot(ax, fi(~qi), ip_tot(~qi), 'x', 'Color', ink, 'MarkerSize', 10, 'LineWidth', 1.5, 'HandleVisibility', 'off'); end
if any(~qc), plot(ax, fc(~qc), ipc_main(~qc), 'x', 'Color', ink, 'MarkerSize', 10, 'LineWidth', 1.5, 'HandleVisibility', 'off'); end
if ~isempty(ph.free), text(ax, mean(ph.free), ymax - 0.7, 'free flight', 'Color', ink2, 'HorizontalAlignment', 'center', 'FontSize', 10); end
if ~isempty(ph.impact), text(ax, mean(ph.impact), ymax - 0.7, 'impact', 'Color', ink2, 'HorizontalAlignment', 'center', 'FontSize', 10); end
if ~isempty(ph.settling), text(ax, mean(ph.settling), ymax - 0.7, 'settling', 'Color', ink2, 'HorizontalAlignment', 'center', 'FontSize', 10); end
if ~isempty(Dsel) && Dsel.frame <= numel(ip_tot)
    plot(ax, Dsel.frame, ip_tot(Dsel.frame) + 0.35, 'v', 'Color', ink, 'MarkerFaceColor', ink, 'MarkerSize', 6, 'HandleVisibility', 'off');
    text(ax, Dsel.frame - 1.5, ip_tot(Dsel.frame) + 0.5, sprintf('replayed frame %d', Dsel.frame), 'Color', ink2, ...
        'HorizontalAlignment', 'right', 'FontSize', 9);
end
xlim(ax, [0.5, nf + 0.5]); ylim(ax, [0, ymax]);
style(ax, ink, ink2);
ylabel(ax, 'full-size solves in the frame');
legend(ax, {sprintf('hard IP main only (qualified total %d)%s', sum(ip_main(qi)), tag(Rip)), ...
    sprintf('IPC main systems (qualified total %d)%s', sum(ipc_main(qc)), tag(Ripc)), ...
    sprintf('hard IP main + verification (qualified total %d)%s', sum(ip_tot(qi)), tag(Rip))}, ...
    'Location', 'northeast', 'Box', 'off', 'TextColor', ink, 'Interpreter', 'none');
title(ax, 'Newton solves per frame, each solver on its own trajectory', 'Color', ink, 'FontWeight', 'normal', 'FontSize', 12);
ax2 = nexttile(tl1); hold(ax2, 'on');
cmax = max([sum(ip_tot), sum(ipc_main), 1]) * 1.08;
if ~isempty(ph.impact), draw_band(ax2, ph.impact, cmax, band); end
plot(ax2, fi, cumsum(ip_main), '--', 'Color', cIP, 'LineWidth', 1);
plot(ax2, fc, cumsum(ipc_main), '-', 'Color', cIPC, 'LineWidth', 2);
plot(ax2, fi, cumsum(ip_tot), '-', 'Color', cIP, 'LineWidth', 2);
text(ax2, fi(end) - 1, sum(ip_tot), sprintf('IP main + verification  %d%s', sum(ip_tot), tag(Rip)), 'Color', ink, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', 'FontSize', 9, 'Interpreter', 'none');
text(ax2, fi(end) - 1, sum(ip_main), sprintf('IP main  %d', sum(ip_main)), 'Color', ink, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', 'FontSize', 9);
text(ax2, fc(end) - 1, sum(ipc_main), sprintf('IPC  %d%s', sum(ipc_main), tag(Ripc)), 'Color', ink, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 9, 'Interpreter', 'none');
xlim(ax2, [0.5, nf + 0.5]); ylim(ax2, [0, cmax]);
style(ax2, ink, ink2);
xlabel(ax2, 'frame'); ylabel(ax2, 'cumulative full-size solves (attempted frames)');
files{end+1} = save_fig(f1, fullfile(outdir, 'newton_solves_vs_frame'));

% ------------------------------------------------------------------ 2. selected impact frame (and a supplementary one)
if ~isempty(Dsel), files{end+1} = frame_plot(Dsel, outdir, cIP, cIPC, ink, ink2); end
if nargin >= 6 && ~isempty(Dextra), files{end+1} = frame_plot(Dextra, outdir, cIP, cIPC, ink, ink2); end
close(f1);
end
function file = frame_plot(Dsel, outdir, cIP, cIPC, ink, ink2)
ii = Dsel.ip_iters; ic = Dsel.ipc_iters;
im = ii(strcmp({ii.kind}, 'main')); iv = ii(strcmp({ii.kind}, 'verification'));
n_ip = numel(im);
sol_ip = 1:n_ip; sol_ver = n_ip + (1:numel(iv));
sol_ipc = 1:numel(ic);
f2 = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1100 640]);
tl = tiledlayout(f2, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
a1 = nexttile(tl); hold(a1, 'on');
semilogy(a1, sol_ipc, [ic.rd_inf], '-o', 'Color', cIPC, 'MarkerFaceColor', cIPC, 'LineWidth', 2, 'MarkerSize', 7);
semilogy(a1, sol_ip, [im.rd_inf], '-s', 'Color', cIP, 'MarkerFaceColor', cIP, 'LineWidth', 2, 'MarkerSize', 7);
if ~isempty(iv)
    semilogy(a1, sol_ver, [iv.rd_inf], 's', 'Color', cIP, 'MarkerFaceColor', 'w', 'LineWidth', 1.5, 'MarkerSize', 8);
end
yline(a1, 5, '--', 'target 5 F_w', 'Color', ink2, 'LabelHorizontalAlignment', 'left', 'FontSize', 9, 'HandleVisibility', 'off');
set(a1, 'YScale', 'log');
style(a1, ink, ink2);
ylabel(a1, '||r_d||_\infty at the iterate (F_w)');
lg = {sprintf('IPC at IP''s frame-start state (%d main)', numel(ic)), sprintf('hard IP main (%d)', n_ip)};
if ~isempty(iv), lg{end+1} = sprintf('hard IP verification (%d)', numel(iv)); end
legend(a1, lg, 'Location', 'northeast', 'Box', 'off', 'TextColor', ink);
title(a1, sprintf('Frame %d, both solvers replayed from IP''s frame-start state (IP keeps its multiplier history)', Dsel.frame), ...
    'Color', ink, 'FontWeight', 'normal', 'FontSize', 12, 'Interpreter', 'none');
xlim(a1, [0.5, max([sol_ver, sol_ip, sol_ipc]) + 0.5]);
a2 = nexttile(tl); hold(a2, 'on');
aipc = [ic.alpha_eff]; aipc(strcmp({ic.limiter}, 'stop')) = NaN;
aip = [im.alpha_eff];
plot(a2, sol_ipc, aipc, '-o', 'Color', cIPC, 'MarkerFaceColor', cIPC, 'LineWidth', 2, 'MarkerSize', 7);
plot(a2, sol_ip, aip, '-s', 'Color', cIP, 'MarkerFaceColor', cIP, 'LineWidth', 2, 'MarkerSize', 7);
% label every shortened IP step with its limiter (s: slack, m: multiplier, e: element, b: merit backtracking,
% g: ground fallback) and, when a contact row set it, the vertex
for q = 1:n_ip
    if aip(q) < 1 && ~isempty(im(q).limiter)
        lb = limiter_tag(im(q).limiter);
        if isfinite(im(q).lim_vertex), lb = sprintf('%s v%d', lb, im(q).lim_vertex); end
        text(a2, sol_ip(q), aip(q) - 0.06, lb, 'Color', ink2, 'FontSize', 8, 'HorizontalAlignment', 'center');
    end
end
ylim(a2, [0, 1.08]); xlim(a2, a1.XLim);
style(a2, ink, ink2);
xlabel(a2, 'solve number within the frame'); ylabel(a2, 'accepted position step fraction');
legend(a2, {'IPC at IP''s state (last system: matched stop, no step)', 'hard IP (labels: limiter, limiting vertex)'}, 'Location', 'southeast', ...
    'Box', 'off', 'TextColor', ink);
file = save_fig(f2, fullfile(outdir, sprintf('impact_frame_%d_residual_and_step', Dsel.frame)));
close(f2);
end

function draw_band(ax, frames, ymax, c)
patch(ax, [frames(1)-0.5 frames(end)+0.5 frames(end)+0.5 frames(1)-0.5], [0 0 ymax ymax], c, ...
    'EdgeColor', 'none', 'HandleVisibility', 'off');
end

function style(ax, ink, ink2)
set(ax, 'Box', 'off', 'TickDir', 'out', 'XColor', ink2, 'YColor', ink2, 'FontSize', 10, 'LineWidth', 0.75, ...
    'GridColor', [0.85 0.85 0.85], 'GridAlpha', 1, 'Layer', 'top');
grid(ax, 'on'); ax.XGrid = 'off';
ax.XLabel.Color = ink; ax.YLabel.Color = ink;
end

function f = save_fig(fig, base)
exportgraphics(fig, [base '.png'], 'Resolution', 150);
savefig(fig, [base '.fig']);
f = [base '.png'];
end

function c = hex(h)
c = sscanf(h(2:end), '%2x%2x%2x', [1 3]) / 255;
end

function v = ternary(c, a, b)
if c, v = a; else, v = b; end
end

function t = limiter_tag(l)
if contains(l, 'slack'), t = 's'; elseif contains(l, 'multiplier'), t = 'm'; elseif contains(l, 'element'), t = 'e';
elseif contains(l, 'merit'), t = 'b'; elseif contains(l, 'rejected'), t = 'r'; else, t = ''; end
if contains(l, 'ground fallback'), t = [t '+g']; end
end
