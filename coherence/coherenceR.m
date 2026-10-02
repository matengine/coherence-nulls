function R = coherenceR(X, Y)
%COHERENCER Unsquared magnitude coherence between two spectral arrays.
%   R = COHERENCER(X, Y) returns the magnitude of the coherency, 1 x nFreq.
%
%   This is R, not R squared. Reported as one or the other under the same name,
%   a value of 0.25 means 6% or 25% of linearly predictable power, and the
%   convention usually follows the analysis toolbox rather than a choice.

R = abs(coherence.coherency(X, Y));
end
