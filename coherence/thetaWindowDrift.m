function [sdHz, meanHz, nWindows] = thetaWindowDrift(x, fs, band, winS, keep)
%THETAWINDOWDRIFT Variability of a slow rhythm's frequency between windows.
%   [SDHZ, MEANHZ, N] = THETAWINDOWDRIFT(X, FS) band-passes x over 5 to 11 Hz,
%   takes the derivative of the unwrapped Hilbert phase as the instantaneous
%   frequency, averages it within one-second windows, and returns the standard
%   deviation of those window means.
%   THETAWINDOWDRIFT(X, FS, BAND, WINS, KEEP) sets the band, the window length
%   and a logical mask over windows, for example fast-running windows only.
%
%   The window mean is used rather than the raw instantaneous frequency because
%   what smears a comb in a pooled spectrum is the frequency difference between
%   analysis segments, not the phase noise inside one cycle.
%
%   Apply this identically to model output and to recordings, otherwise the two
%   are not on the same axis and the comparison means nothing.

if nargin < 3 || isempty(band), band = [5 11]; end
if nargin < 4 || isempty(winS), winS = 1.0; end
if nargin < 5, keep = []; end
x = x(:);
[b, a] = butter(3, band / (fs / 2), 'bandpass');
ph = unwrap(angle(hilbert(filtfilt(b, a, x - mean(x)))));
inst = diff(ph) * fs / (2 * pi);
n = round(winS * fs);
nw = floor(numel(inst) / n);
fw = mean(reshape(inst(1:nw * n), n, nw), 1).';
m = isfinite(fw) & fw > band(1) - 1 & fw < band(2) + 1;
if ~isempty(keep)
    k = logical(keep(:));
    k = k(1:min(nw, numel(k)));
    k(end + 1:nw) = false;
    m = m & k;
end
if sum(m) < 3
    sdHz = NaN; meanHz = NaN; nWindows = sum(m);
    return
end
sdHz = std(fw(m));
meanHz = mean(fw(m));
nWindows = sum(m);
end
