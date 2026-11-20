function ph = trajectory_phases(Rip, Ripc, P)
%TRAJECTORY_PHASES  Free flight / impact / settling, from the first frame in which either solver has a contact pair.
% Defined over the requested frames 1..P.frames, independent of where a trajectory stopped. Impact = the first
% contact frame and the next P.impact_frames - 1 frames (the reference's 20-frame window); settling = the rest.
nf = P.frames;
c1 = find([Rip.frames.contact], 1); c2 = find([Ripc.frames.contact], 1);
fc = min([c1, c2]);
ph.rule = sprintf(['free flight = frames before the first frame in which either solver has a contact pair; impact = ' ...
    'that frame and the next %d; settling = the rest'], P.impact_frames - 1);
if isempty(fc)
    ph.first_contact = NaN; ph.first_contact_ip = NaN; ph.first_contact_ipc = NaN;
    ph.free = 1:nf; ph.impact = zeros(1, 0); ph.settling = zeros(1, 0); ph.all = 1:nf;
    return;
end
ph.first_contact = fc; ph.first_contact_ip = empty_nan(c1); ph.first_contact_ipc = empty_nan(c2);
ph.free = 1:fc-1;
ph.impact = fc:min(nf, fc + P.impact_frames - 1);
ph.settling = ph.impact(end)+1:nf;
ph.all = 1:nf;
end

function v = empty_nan(v)
if isempty(v), v = NaN; end
end
