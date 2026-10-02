function bad = badChannels(X, fs, band, window, corrMin, varRatio)
%BADCHANNELS Channels that are dead, grossly noisy, or unrelated to neighbours.
%   BAD = BADCHANNELS(X, FS) returns the indices of channels to interpolate
%   before a spatial derivative.
%
%   The criterion is permissive on purpose. A strong local generator
%   legitimately decorrelates a channel from its neighbours, so only channels
%   essentially unrelated to both neighbours, or with variance an order of
%   magnitude away from the median, are flagged. Over-flagging destroys real
%   laminar structure, which is worse than leaving a marginal channel in.

if nargin < 3 || isempty(band), band = [1 100]; end
if nargin < 4 || isempty(window), window = [600 1200]; end
if nargin < 5 || isempty(corrMin), corrMin = 0.35; end
if nargin < 6 || isempty(varRatio), varRatio = 10; end
X = double(X);
seg = recordings.analysisSlice(size(X, 2), fs, window);
[b, a] = butter(2, band / (fs / 2), 'bandpass');
Xb = filtfilt(b, a, X(:, seg).').';
C = corrcoef(Xb.');
n = size(X, 1);
v = var(Xb, 1, 2);
mv = median(v);
bad = [];
for i = 1:n
    nb = [i - 1, i + 1];
    nb = nb(nb >= 1 & nb <= n);
    if max(C(i, nb)) < corrMin || v(i) > varRatio * mv || v(i) < mv / varRatio
        bad(end + 1) = i; %#ok<AGROW>
    end
end
end
