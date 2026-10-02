function out = troughOffsets(xa, xb, fs, band, keepSamples, maxAbsMs)
%TROUGHOFFSETS Cross-regional timing of a slow rhythm, trough to nearest trough.
%   OUT = TROUGHOFFSETS(XA, XB, FS) returns the intervals in milliseconds from
%   each trough of xa to the nearest trough of xb, with both summaries.
%
%   Two summaries are returned because they answer different questions and can
%   differ when the distribution has heavy tails: the linear standard deviation
%   of the intervals, and the circular standard deviation from the resultant
%   length at the band centre.
%
%   This is the measurement that decides whether a shared pacemaker could have
%   aligned two independent fast generators, which needs sub-millisecond
%   precision.

if nargin < 4 || isempty(band), band = [5 11]; end
if nargin < 5, keepSamples = []; end
if nargin < 6 || isempty(maxAbsMs), maxAbsMs = 62; end
[b, a] = butter(2, band / (fs / 2), 'bandpass');
ta = troughsOf(xa, b, a);
tb = troughsOf(xb, b, a);
if ~isempty(keepSamples)
    k = logical(keepSamples(:));
    ta = ta(k(ta));
end
if numel(tb) < 2 || isempty(ta)
    out = struct('dtMs', [], 'n', 0, 'meanMs', NaN, 'sdMs', NaN, 'resultant', NaN, 'circularSdDeg', NaN);
    return
end
near = tb(interp1(tb, 1:numel(tb), ta, 'nearest', 'extrap'));
dt = (near - ta) / fs * 1000;
dt = dt(abs(dt) < maxAbsMs);
fMid = 0.5 * (band(1) + band(2));
Rbar = abs(mean(exp(1i * 2 * pi * fMid * dt / 1000)));
out = struct('dtMs', dt, 'n', numel(dt), 'meanMs', mean(dt), 'sdMs', std(dt, 1), ...
             'resultant', Rbar, 'circularSdDeg', coherence.sigmaFromR(Rbar));
end

function idx = troughsOf(x, b, a)
p = angle(hilbert(filtfilt(b, a, double(x(:)))));
w = coherence.wrapAngle(p - pi);
idx = find(w(1:end - 1) < 0 & w(2:end) >= 0) + 1;
end
