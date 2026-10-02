function [Leff, gA, gB] = lEffWithinBin(X, Y, index, nFreq, maxLag)
%LEFFWITHINBIN Effective estimate count when only a subset of segments is used.
%   [LEFF, GA, GB] = LEFFWITHINBIN(X, Y, INDEX, NFREQ) measures persistence
%   using only pairs whose partner at lag j is also retained, so it is measured
%   on segments genuinely consecutive in the recording rather than across the
%   gaps that subsetting creates.
%
%   Measuring persistence on a non-contiguous subset is a real trap: it reports
%   almost no persistence where there is plenty.

if nargin < 5, maxLag = 10; end
index = index(:);
L = numel(index);
gA = laggedSubset(X, index, nFreq, maxLag);
gB = laggedSubset(Y, index, nFreq, maxLag);
j = (1:maxLag).';
s = sum((1 - j / L) .* real(gA .* conj(gB)), 1, 'omitnan');
Leff = L ./ (1 + 2 * s);
end

function g = laggedSubset(Z, index, nFreq, maxLag)
g = nan(maxLag, nFreq);
den = mean(abs(Z(index, :)) .^ 2, 1);
for j = 1:maxLag
    k = index(ismember(index + j, index));
    if numel(k) < 20
        continue
    end
    g(j, :) = mean(Z(k, :) .* conj(Z(k + j, :)), 1) ./ den;
end
end
