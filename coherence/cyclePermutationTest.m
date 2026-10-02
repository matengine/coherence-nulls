function [obs, sur, f, Robs] = cyclePermutationTest(xa, xb, troughs, fs, band, nSurr, nfft, stat)
%CYCLEPERMUTATIONTEST Conditional null: more than shared rhythmic drive?
%   [OBS, SUR, F, ROBS] = CYCLEPERMUTATIONTEST(XA, XB, TROUGHS, FS, BAND, NSURR)
%   returns the observed same-cycle band statistic and its surrogate
%   distribution. Compare OBS with coherence.quantileLinear(SUR, 0.95).
%
%   Each cycle of one region is paired with a different cycle of the other,
%   keeping its internal structure and its position relative to the trough. A
%   high surrogate is expected: locking of both regions to the shared rhythm
%   survives permutation.
%
%   The surrogate depends on the boundaries sitting on the event that resets
%   the fast rhythm. Boundary error only lowers it, by exp(-(2*pi*f*sigma)^2),
%   which at 80 Hz is 0.78 at 1 ms and 0.02 at 4 ms, so the failure is
%   one-sided: a non-exceedance is valid with any boundaries, an exceedance
%   only with boundaries good to about a millisecond. Boundaries from a
%   receiver's own field do not qualify when its slow rhythm mixes shared
%   drive with local components. Calibrate on a shared-drive simulation first.
%
%   Working on epochs rather than a reassembled trace avoids two artefacts:
%   time-warping onto a uniform phase grid slows a burst in a stretched cycle
%   and manufactures false positives, and splicing cycles of differing length
%   injects broadband edge energy that lowers the surrogate.

if nargin < 7, nfft = []; end
if nargin < 8 || isempty(stat), stat = 'mean'; end
[f, EA, EB, durs] = coherence.cycleEpochFFTs(xa, xb, troughs, fs, nfft);
Robs = coherence.coherenceR(EA, EB);
obs = pick(Robs, f, band, stat);
sur = zeros(nSurr, 1);
for k = 1:nSurr
    Rs = coherence.coherenceR(EA, coherence.permuteEpochsByDuration(EB, durs));
    sur(k) = pick(Rs, f, band, stat);
end
end

function v = pick(R, f, band, stat)
if strcmpi(stat, 'peak')
    v = coherence.bandPeak(R, f, band);
else
    v = coherence.bandMean(R, f, band);
end
end
