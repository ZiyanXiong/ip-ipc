function [first, nonpos_start, dets] = elem_check(S, x, dx, amax)
%ELEM_CHECK  Shared element-validity check of the step x + t dx, t in (0, amax].
%   first        smallest t in (0, amax] at which some triangle's det(Ds) reaches zero (Inf if none)
%   nonpos_start number of triangles with det(Ds) <= 0 at t = 0
%   dets         det(Ds) at t = 0 (per triangle)
% det(Ds + t dDs) = c0 + c1 t + c2 t^2 is exact for linear triangles.
T = S.T;
Ds11 = x(2*T(:,2)-1) - x(2*T(:,1)-1); Ds21 = x(2*T(:,2)) - x(2*T(:,1));
Ds12 = x(2*T(:,3)-1) - x(2*T(:,1)-1); Ds22 = x(2*T(:,3)) - x(2*T(:,1));
c0 = Ds11 .* Ds22 - Ds12 .* Ds21;
dets = c0;
nonpos_start = sum(~(c0 > 0));
first = Inf;
if nargin < 3 || isempty(dx), return; end
if nargin < 4, amax = 1; end
d11 = dx(2*T(:,2)-1) - dx(2*T(:,1)-1); d21 = dx(2*T(:,2)) - dx(2*T(:,1));
d12 = dx(2*T(:,3)-1) - dx(2*T(:,1)-1); d22 = dx(2*T(:,3)) - dx(2*T(:,1));
c2 = d11 .* d22 - d12 .* d21;
c1 = Ds11 .* d22 + d11 .* Ds22 - Ds12 .* d21 - d12 .* Ds21;
t = inf(size(c0));
scale = abs(c0) + abs(c1) + abs(c2);
lin = abs(c2) <= 1e-14 * scale;
il = lin & c1 < 0;
t(il) = -c0(il) ./ c1(il);
q = ~lin;
disc = c1.^2 - 4 * c2 .* c0;
iq = q & disc >= 0;
sq = sqrt(max(disc(iq), 0));
qq = -0.5 * (c1(iq) + sign(c1(iq) + (c1(iq) == 0)) .* sq);
r1 = qq ./ c2(iq);
r2 = c0(iq) ./ qq;
r1(~(r1 > 0)) = Inf; r2(~(r2 > 0)) = Inf;
t(iq) = min(r1, r2);
t(~(t > 0)) = Inf;
first = min(t);
if first > amax, first = Inf; end
end
