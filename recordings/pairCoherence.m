function out = pairCoherence(Afull, Bfull, keep, f, nSurr, bands)
%PAIRCOHERENCE Coherence of one channel pair in one speed bin, with its nulls.
%   OUT = PAIRCOHERENCE(AFULL, BFULL, KEEP, F) where AFULL and BFULL are
%   whole-session (segment, frequency) arrays and KEEP selects this bin.
%   OUT = PAIRCOHERENCE(..., NSURR, BANDS) sets the surrogate count and the
%   bands, a cell array of {name, [lo hi], 'mean'|'peak'}.
%
%   The circular shift is applied to B's whole sequence before the mask;
%   shifting a subset would reshuffle data that are not contiguous in time.
%
%   Returns R, the per-frequency 95th percentiles of both nulls, L, Leff, the
%   lagged self-coherences, and each band statistic with its own surrogate
%   quantile, bias-corrected effect size and equivalent phase scatter.

if nargin < 5 || isempty(nSurr), nSurr = 200; end
if nargin < 6 || isempty(bands)
    bands = {'theta', [5 11], 'peak'; 'gamma', [60 100], 'mean'};
end
keep = logical(keep(:));
index = find(keep);
L = numel(index);
N = size(Bfull, 1);
A = Afull(keep, :);
R = coherence.coherenceR(A, Bfull(keep, :));

shifts = randi([5, N - 5], 1, nSurr);
Rshift = zeros(nSurr, size(A, 2));
for k = 1:nSurr
    Ys = circshift(Bfull, shifts(k), 1);
    Rshift(k, :) = coherence.coherenceR(A, Ys(keep, :));
end
Rperm = coherence.nullPermutation(A, Bfull(keep, :), max(round(nSurr / 2), 10));

out = struct();
out.R = R;
out.q95shift = coherence.quantileLinear(Rshift, 0.95);
out.q95perm = coherence.quantileLinear(Rperm, 0.95);
out.L = L;
[out.Leff, out.laggedA, out.laggedB] = coherence.lEffWithinBin(Afull, Bfull, index, numel(f));

for i = 1:size(bands, 1)
    name = bands{i, 1};
    band = bands{i, 2};
    stat = bands{i, 3};
    if strcmpi(stat, 'peak')
        obs = coherence.bandPeak(R, f, band);
        sur = coherence.bandPeak(Rshift, f, band);
    else
        obs = coherence.bandMean(R, f, band);
        sur = coherence.bandMean(Rshift, f, band);
    end
    m = coherence.bandMask(f, band);
    fb = f(m);
    Rb = R(m);
    [~, j] = max(Rb);
    [~, iBand] = min(abs(f - fb(j)));
    Lhere = max(out.Leff(iBand), 2);
    r2 = coherence.biasCorrectedR2(obs, Lhere);
    s = struct('statistic', stat, 'observed', obs, ...
               'null95', coherence.quantileLinear(sur, 0.95), 'nullMean', mean(sur), ...
               'Leff', Lhere, 'R2corrected', r2);
    if obs > 0
        s.sigmaDeg = coherence.sigmaFromR(sqrt(r2));
    else
        s.sigmaDeg = Inf;
    end
    out.bands.(name) = s;
end
end
