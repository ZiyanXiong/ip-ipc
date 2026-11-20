function [B, dB, ddB] = ipc_barrier(D, Dhat)
%IPC_BARRIER  The native IPC barrier in the squared-distance convention, without kappa (zero thickness).
%   codim_ipc_contact.inl  KappaBarrierLog2 / dKappaBarrierLog2dD / ddKappaBarrierLog2ddD  (xi == 0 branch):
%   B(D)    = (D - Dhat)^2 log(D/Dhat)^2,                      D = d^2, Dhat = d_hat^2, active for D < Dhat
%   dB/dD   = (2D - 2Dhat) log^2 + 2 (D - Dhat)^2 log / D
%   d2B/dD2 = 2 log^2 - log * 2 (D-Dhat)^2/D^2 + 2 (D-Dhat)^2/D^2 + 4 log (2D - 2Dhat) / D
% Outside the band (D >= Dhat) all three are zero.
B = zeros(size(D)); dB = B; ddB = B;
a = D < Dhat;
Da = D(a);
x1 = log(Da / Dhat);
B(a) = (Da - Dhat).^2 .* x1.^2;
dB(a) = x1.^2 .* (2*Da - 2*Dhat) + 2 * x1 .* (Da - Dhat).^2 ./ Da;
x3 = 2 * (Da - Dhat).^2 ./ Da.^2;
ddB(a) = 2 * x1.^2 - x1 .* x3 + x3 + 4 * x1 .* (2*Da - 2*Dhat) ./ Da;
end
