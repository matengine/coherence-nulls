function [f, F, speedSeg, ok] = segmentSession(channels, fs, speedPerSample, segLenS, artifactSd)
%SEGMENTSESSION Cut named channels into segments and return their spectra.
%   [F, FSPEC, SPEEDSEG, OK] = SEGMENTSESSION(CHANNELS, FS, SPEED) takes a
%   struct whose fields are 1-D traces of equal length and returns a struct of
%   (segment, frequency) arrays with the same field names, the mean running
%   speed of each segment, and a mask of segments free of artifact and with
%   tracking.
%
%   Segments are one second, non-overlapping, demeaned and Hann-tapered.

if nargin < 4 || isempty(segLenS), segLenS = 1.0; end
if nargin < 5 || isempty(artifactSd), artifactSd = 10; end
names = fieldnames(channels);
n = round(segLenS * fs);
nsamp = inf;
for i = 1:numel(names)
    nsamp = min(nsamp, numel(channels.(names{i})));
end
N = floor(nsamp / n);
sp = speedPerSample(:);
speedSeg = mean(reshape(sp(1:N * n), n, N), 1, 'omitnan').';
ok = ~isnan(speedSeg);
for i = 1:numel(names)
    x = channels.(names{i});
    xs = reshape(x(1:N * n), n, N).';
    ok = ok & (max(abs(xs - mean(xs(:))), [], 2) < artifactSd * std(xs(:), 1));
end
F = struct();
for i = 1:numel(names)
    x = channels.(names{i});
    [f, F.(names{i})] = coherence.segmentFFTs(x(1:N * n), fs, segLenS);
end
end
