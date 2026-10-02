function d = wrappedNormalPdf(phiDeg, sigmaDeg, muDeg, nWrap)
%WRAPPEDNORMALPDF Density per degree of a wrapped normal.

if nargin < 3, muDeg = 0; end
if nargin < 4, nWrap = 4; end
phi = deg2rad(phiDeg);
s = deg2rad(sigmaDeg);
mu = deg2rad(muDeg);
d = zeros(size(phi));
for k = -nWrap:nWrap
    d = d + exp(-0.5 * ((phi - mu + 2 * pi * k) / s) .^ 2);
end
d = d / (s * sqrt(2 * pi)) * pi / 180;
end
