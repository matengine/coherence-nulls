function Rs = nullPermutation(X, Y, nSurr)
%NULLPERMUTATION Detection null by random re-pairing of segments.
%   Kept for comparison. This destroys each signal's ordering as well as the
%   pairing, so the surrogate spectra lose the serial persistence the real data
%   have. The resulting null is too low wherever a rhythm persists and the test
%   is anti-conservative there. Prefer coherence.nullCircularShift.

N = size(X, 1);
Rs = zeros(nSurr, size(X, 2));
for k = 1:nSurr
    Rs(k, :) = coherence.coherenceR(X, Y(coherence.derangement(N), :));
end
end
