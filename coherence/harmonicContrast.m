function h = harmonicContrast(Rf, f, band, L, fSlow, halfwidth)
%HARMONICCONTRAST Coherence on the harmonics minus coherence between them.
%   H = HARMONICCONTRAST(RF, F, BAND, L, FSLOW) takes the mean of R within
%   1 Hz of a multiple of FSLOW, subtracts the mean within 1 Hz of a
%   half-integer multiple, and divides by the sampling standard deviation.
%   H = HARMONICCONTRAST(..., HALFWIDTH) sets that window.
%
%   Large and positive when a slow rhythm periodically resets the fast
%   generator, because a reset is cyclostationary at the slow period and places
%   its coherence on the harmonics. About zero for a band-limited generator,
%   whatever its comb index. NaN if either set has fewer than two bins.

if nargin < 6, halfwidth = 1.0; end
f = f(:).';
Rf = Rf(:).';
m = coherence.bandMask(f, band);
on = false(size(f));
off = false(size(f));
for k = 1:(floor(max(f) / max(fSlow, 1e-6)) + 1)
    on = on | abs(f - k * fSlow) <= halfwidth;
    off = off | abs(f - (k + 0.5) * fSlow) <= halfwidth;
end
on = on & m;
off = off & m;
if sum(on) < 2 || sum(off) < 2
    h = NaN;
    return
end
r = Rf(m);
h = (mean(Rf(on)) - mean(Rf(off))) / ((1 - mean(r) ^ 2) / sqrt(2 * L));
end
