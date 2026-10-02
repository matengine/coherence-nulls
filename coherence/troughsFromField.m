function idx = troughsFromField(x, fs, band)
%TROUGHSFROMFIELD Cycle boundaries from a band-passed field potential.
%   This is the ordinary procedure, and the one that makes the conditional null
%   unusable for a positive claim: when a field's slow rhythm mixes a shared
%   drive with independent local components, these boundaries sit many
%   milliseconds from the event that resets the fast rhythm. Use them for a
%   non-exceedance result, or to demonstrate the failure, not to support a
%   claim about a channel.

if nargin < 3, band = [5 11]; end
[b, a] = butter(2, band / (fs / 2), 'bandpass');
idx = coherence.troughsFromPhase(angle(hilbert(filtfilt(b, a, x(:)))));
end
