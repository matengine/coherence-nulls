function out = currentSourceDensity(X, bad)
%CURRENTSOURCEDENSITY Negative second spatial difference along the channel axis.
%   OUT = CURRENTSOURCEDENSITY(X) with X as (channel, sample).
%   OUT = CURRENTSOURCEDENSITY(X, BAD) interpolates the listed channel indices
%   from their nearest good neighbours first.
%
%   Sink-negative convention. The first and last channels have no CSD and come
%   back as NaN. Assumes channel order follows depth with uniform spacing,
%   which the correlation matrix from recordings.laminarProfile will show.
%
%   This step decides what a cross-regional coherence means. A field-potential
%   coherence with no spectral peak that rises toward the Nyquist band is a
%   shared broadband component, not a channel between regions, and it can
%   reach 0.9 where the CSD at the same sites sits near zero.

if nargin < 2, bad = []; end
Xi = double(X);
n = size(Xi, 1);
bad = bad(:).';
for i = bad
    lo = i - 1;
    while lo >= 1 && any(bad == lo), lo = lo - 1; end
    hi = i + 1;
    while hi <= n && any(bad == hi), hi = hi + 1; end
    if lo >= 1 && hi <= n
        Xi(i, :) = 0.5 * (Xi(lo, :) + Xi(hi, :));
    elseif lo >= 1
        Xi(i, :) = Xi(lo, :);
    elseif hi <= n
        Xi(i, :) = Xi(hi, :);
    end
end
out = nan(size(Xi));
out(2:end - 1, :) = -(Xi(1:end - 2, :) - 2 * Xi(2:end - 1, :) + Xi(3:end, :));
end
