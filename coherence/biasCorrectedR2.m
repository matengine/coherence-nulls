function r2 = biasCorrectedR2(R, L)
%BIASCORRECTEDR2 Approximately unbiased population R^2, truncated at zero.
%   The raw estimate is biased upward by roughly 1/L, which is large at the
%   small effective sample sizes that serial persistence produces at a rhythm's
%   own frequency.

r2 = max((R .^ 2 - 1 ./ L) ./ (1 - 1 ./ L), 0);
end
