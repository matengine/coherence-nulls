function [Leff, gxx, gyy] = lEffFrequency(X, Y, maxLag)
%LEFFFREQUENCY Independent spectral estimates per frequency, not segment count.
%   [LEFF, GXX, GYY] = LEFFFREQUENCY(X, Y, MAXLAG) returns
%   N / (1 + 2*sum_j (1 - j/N) Re[gxx(f,j) conj(gyy(f,j))]).
%
%   Segments are not independent at a frequency where a rhythm persists across
%   them, so pass this, not the segment count, to coherence.critR,
%   coherence.pValue and coherence.biasCorrectedR2.
%
%   Cap it at N before using it for inference: a value above N is a variance
%   ratio in a finite sample, not a count of independent observations.

if nargin < 3, maxLag = 30; end
N = size(X, 1);
gxx = coherence.laggedSelfCoherence(X, maxLag);
gyy = coherence.laggedSelfCoherence(Y, maxLag);
j = (1:maxLag).';
s = sum((1 - j / N) .* real(gxx .* conj(gyy)), 1);
Leff = N ./ (1 + 2 * s);
end
