function [q, T] = bandNull(Rsurr, f, band, p, stat)
%BANDNULL Null distribution of a band statistic, not of a single frequency.
%   [Q, T] = BANDNULL(RSURR, F, BAND, P, STAT) returns the P quantile of the
%   band statistic across surrogates and the full surrogate distribution.
%
%   STAT must match whatever is reported: 'mean' for a band average, 'peak' for
%   a band maximum. Averaging the bins of a band shrinks its null sharply, by
%   roughly the square root of the number of independent bins, so a
%   per-frequency threshold is not a substitute.

if nargin < 4 || isempty(p), p = 0.95; end
if nargin < 5 || isempty(stat), stat = 'mean'; end
if strcmpi(stat, 'peak')
    T = coherence.bandPeak(Rsurr, f, band);
else
    T = coherence.bandMean(Rsurr, f, band);
end
q = coherence.quantileLinear(T, p);
end
