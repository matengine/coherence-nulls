function [f, EA, EB, durs] = cycleEpochFFTs(xa, xb, troughs, fs, nfft, minLen, maxLen)
%CYCLEEPOCHFFTS One spectral estimate per cycle, from epochs between troughs.
%   [F, EA, EB, DURS] = CYCLEEPOCHFFTS(XA, XB, TROUGHS, FS) demeans each cycle,
%   Hann-tapers it over its own length and zero-pads to a common transform
%   length. EA and EB are nCycles x nFreq.
%
%   The taper puts low weight at the epoch edges, so a shared event sitting on
%   a boundary is down-weighted in the observed statistic. That changes the
%   magnitude but not the verdict; a rectangular window raises the observation
%   and its surrogate together.

if nargin < 5 || isempty(nfft), nfft = round(fs); end
if nargin < 6, minLen = []; end
if nargin < 7, maxLen = []; end
troughs = troughs(:);
durs = diff(troughs);
keep = true(size(durs));
if ~isempty(minLen), keep = keep & durs >= minLen; end
if ~isempty(maxLen), keep = keep & durs <= maxLen; end
starts = troughs(1:end - 1);
starts = starts(keep);
durs = durs(keep);
nf = floor(nfft / 2) + 1;
EA = zeros(numel(durs), nf);
EB = zeros(numel(durs), nf);
for j = 1:numel(durs)
    d = durs(j);
    s0 = starts(j);
    w = coherence.hannPeriodic(d);
    a = xa(s0:s0 + d - 1); a = a(:) - mean(a);
    b = xb(s0:s0 + d - 1); b = b(:) - mean(b);
    FA = fft(a .* w, nfft);
    FB = fft(b .* w, nfft);
    EA(j, :) = FA(1:nf).';
    EB(j, :) = FB(1:nf).';
end
f = (0:floor(nfft / 2)) * fs / nfft;
end
