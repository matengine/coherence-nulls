function r = critR(L, alpha)
%CRITR Pointwise critical coherence at level alpha for L independent estimates.
%   Follows from the Goodman (1957) null, P(Rhat > r) = (1 - r^2)^(L-1).
%
%   L is the number of INDEPENDENT estimates. Where a rhythm persists across
%   segments that is not the segment count; use coherence.lEffFrequency.

if nargin < 2, alpha = 0.05; end
r = sqrt(1 - alpha .^ (1 ./ (L - 1)));
end
