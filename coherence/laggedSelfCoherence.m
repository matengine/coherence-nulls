function g = laggedSelfCoherence(X, maxLag)
%LAGGEDSELFCOHERENCE How much of a segment's spectrum the next segments repeat.
%   G = LAGGEDSELFCOHERENCE(X, MAXLAG) returns a maxLag x nFreq complex array,
%   <X_k conj(X_{k+j})> / <|X_k|^2> for j = 1..maxLag.
%
%   Its magnitude at lag one is the fraction of a segment's spectral estimate
%   that the following segment repeats, which is what makes segments
%   non-independent at a frequency where a rhythm persists.

N = size(X, 1);
den = mean(abs(X) .^ 2, 1);
g = zeros(maxLag, size(X, 2));
for j = 1:maxLag
    g(j, :) = mean(X(1:N - j, :) .* conj(X(j + 1:N, :)), 1) ./ den;
end
end
