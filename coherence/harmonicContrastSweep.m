function [bestF, bestV] = harmonicContrastSweep(Rf, f, band, L, fLo, fHi, step, halfwidth)
%HARMONICCONTRASTSWEEP Largest harmonic contrast over a sweep of the fundamental.
%   [BESTF, BESTV] = HARMONICCONTRASTSWEEP(RF, F, BAND, L) sweeps the assumed
%   fundamental over 7.5 to 10 Hz in 0.1 Hz steps.
%
%   The statistic depends on which frequency is called the fundamental, and in
%   a recording the coherence peak bin and the mean instantaneous frequency
%   need not agree. The sweep maximum is the value most favourable to a comb,
%   so a small maximum is a strong statement that there is no comb.

if nargin < 5, fLo = 7.5; end
if nargin < 6, fHi = 10.0; end
if nargin < 7, step = 0.1; end
if nargin < 8, halfwidth = 1.0; end
bestF = NaN;
bestV = -Inf;
for t = fLo:step:fHi
    v = coherence.harmonicContrast(Rf, f, band, L, t, halfwidth);
    if isfinite(v) && v > bestV
        bestV = v;
        bestF = t;
    end
end
if ~isfinite(bestV)
    bestV = NaN;
end
end
