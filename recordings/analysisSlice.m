function idx = analysisSlice(nsamp, fs, window, minS)
%ANALYSISSLICE Sample indices for a time window, falling back to the whole run.
%   The default windows elsewhere in this package skip the start of a session,
%   where a recording is often still settling. A short session has no such
%   interior window, so the whole recording is used instead of an empty range,
%   which would otherwise fail deep inside a filter with an unhelpful message.

if nargin < 4, minS = 10; end
lo = max(round(max(window(1), 0) * fs), 1);
hi = min(round(window(2) * fs), nsamp);
if hi - lo < round(minS * fs)
    lo = 1;
    hi = nsamp;
end
idx = lo:hi;
end
