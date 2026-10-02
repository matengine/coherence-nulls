function [f, X] = segmentFFTsMultitaper(x, fs, segLenS, NW, K)
%SEGMENTFFTSMULTITAPER Slepian multitaper spectra, K tapers per segment.
%   [F, X] = SEGMENTFFTSMULTITAPER(X, FS, SEGLENS, NW, K) returns the tapered
%   spectra pooled into (nSegments*K) x nFrequencies.
%
%   Multitaper smoothing averages across neighbouring frequencies. It lowers the
%   variance but does not remove narrow structure the signal actually contains,
%   so a comb survives it in attenuated form rather than disappearing.

if nargin < 3, segLenS = 1.0; end
if nargin < 4, NW = 2.0; end
if nargin < 5, K = 3; end
n = round(segLenS * fs);
N = floor(numel(x) / n);
xs = reshape(x(1:N * n), n, N).';
xs = xs - mean(xs, 2);
tapers = dpss(n, NW, K).';
X = zeros(N * K, floor(n / 2) + 1);
for k = 1:K
    [f, Xk] = coherence.rfftRows(xs .* tapers(k, :), fs);
    X((k - 1) * N + (1:N), :) = Xk;
end
end
