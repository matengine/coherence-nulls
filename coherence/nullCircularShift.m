function Rs = nullCircularShift(X, Y, nSurr, minShift)
%NULLCIRCULARSHIFT Detection null by circular shift of one segment sequence.
%   Each signal keeps its own serial structure and only the correspondence
%   between them is broken, which is what a detection null should do. This is
%   the one to use.

if nargin < 4, minShift = 5; end
N = size(Y, 1);
shifts = randi([minShift, N - minShift], 1, nSurr);
Rs = zeros(nSurr, size(X, 2));
for k = 1:nSurr
    Rs(k, :) = coherence.coherenceR(X, circshift(Y, shifts(k), 1));
end
end
