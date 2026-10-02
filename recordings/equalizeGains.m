function [Xc, gain] = equalizeGains(X, fs, band, window, medianSize)
%EQUALIZEGAINS Per-channel gain correction from the high-frequency noise floor.
%   [XC, GAIN] = EQUALIZEGAINS(X, FS) rescales each channel so that its
%   200 to 400 Hz power follows a running median across depth.
%
%   The noise floor of a laminar probe varies smoothly with depth, so a
%   channel's departure from that trend is treated as a gain mismatch,
%   typically electrode impedance. This matters because a spatial derivative
%   differences neighbouring channels, so an uncorrected gain step becomes a
%   spurious source-sink pair.

if nargin < 3 || isempty(band), band = [200 400]; end
if nargin < 4 || isempty(window), window = [300 1500]; end
if nargin < 5 || isempty(medianSize), medianSize = 7; end
X = double(X);
seg = recordings.analysisSlice(size(X, 2), fs, window);
nper = min(round(fs), numel(seg));
[P, f] = pwelch(X(:, seg).', coherence.hannPeriodic(nper), floor(nper / 2), nper, fs);
sel = f >= band(1) & f <= band(2);
if ~any(sel)
    Xc = X;
    gain = ones(size(X, 1), 1);
    return
end
p = mean(P(sel, :), 1).';
trend = exp(recordings.medianFilterReplicate(log(p), medianSize));
gain = sqrt(trend ./ p);
Xc = X .* gain;
end
