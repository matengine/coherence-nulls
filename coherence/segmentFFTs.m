function [f, X] = segmentFFTs(x, fs, segLenS)
%SEGMENTFFTS Non-overlapping Hann-tapered segments of x, as a spectral array.
%   [F, X] = SEGMENTFFTS(X, FS) cuts x into one-second segments, demeans each,
%   applies a periodic Hann taper and returns the one-sided spectra.
%   [F, X] = SEGMENTFFTS(X, FS, SEGLENS) sets the segment length in seconds.
%
%   X is nSegments x nFrequencies and F is 1 x nFrequencies. The number of rows
%   is the number of spectral estimates entering any later average, which is
%   what makes the sample size explicit everywhere downstream.
%
%   The segment length sets the frequency resolution: 1 s gives 1-Hz bins and
%   0.5 s gives 2-Hz bins. It matters for any band statistic, because the number
%   of bins in a band, and so how far a band average shrinks its own null
%   distribution, follows from it.

if nargin < 3, segLenS = 1.0; end
n = round(segLenS * fs);
N = floor(numel(x) / n);
xs = reshape(x(1:N * n), n, N).';
xs = xs - mean(xs, 2);
w = coherence.hannPeriodic(n).';
[f, X] = coherence.rfftRows(xs .* w, fs);
end
