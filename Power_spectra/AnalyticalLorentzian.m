% AnalyticalLorentzian
% [fc, D, sfc, sD, chi] = AnalyticalLorentzian(spec, df, Tmsr, fstart)
% This function is based on  K. Berg-Sorensen, H. Flyvberg, Rev. Sci.
% Inst., 75, 594 (2004).  It fits a lorentzian to a positive frequency
% power spectra, returning the diffusion constant and the corner frequency
% and the uncertainty in each of these.

% Jeffrey Moffitt
% Corrected chi 2/4/06

function [fc, D, sfc, sD, chi] = AnalyticalLorentzian(spec, df, Tmsr, fstart)

if nargin < 4
    fstart = 0;
end
% Compute useful quantities

f = [0:length(spec) - 1]*df + fstart;

S00 = length(spec);
S01 = sum(spec);
S02 = sum(spec.*spec);
S11 = sum(spec.*f.*f);
S12 = sum(spec.*spec.*f.*f);
S22 = sum(spec.*spec.*f.*f.*f.*f);

fc = sqrt((S01*S22-S11*S12)/(S11*S02 - S01*S12));
D = pi^2*(S02*S22 - S12^2)/(S11*S02 - S01*S12);

chi = S00 - (S01^2*S22 + S11^2*S02 - 2*S01*S11*S12)/(S02*S22-S12^2);

x1 = fstart/fc;
x2 = (fstart + length(spec)*df)/fc;

u = 2*x2/(1+x2^2) -  2*x1/(1+x1^2) + 2*atan((x2-x1)/(1+x2*x1));

v = 4/(x2 - x1)*(atan((x2-x1)/(1+x2*x1)))^2;

sfc = sqrt(fc/((u-v)*Tmsr));

sD = D*sqrt(u/(fc*Tmsr*(u-v)*(x2-x1)));

