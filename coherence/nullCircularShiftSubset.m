function Rs = nullCircularShiftSubset(X, Yfull, keep, nSurr, minShift)
%NULLCIRCULARSHIFTSUBSET Circular-shift null when only a subset is analysed.
%   The shift is applied to the whole sequence and the subset taken afterwards,
%   so the surrogate is built from contiguous data rather than from a
%   reshuffling of an already fragmented subset.

if nargin < 5, minShift = 5; end
N = size(Yfull, 1);
shifts = randi([minShift, N - minShift], 1, nSurr);
Rs = zeros(nSurr, size(X, 2));
for k = 1:nSurr
    Ys = circshift(Yfull, shifts(k), 1);
    Rs(k, :) = coherence.coherenceR(X, Ys(keep, :));
end
end
