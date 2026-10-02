function y = medianFilterReplicate(x, k)
%MEDIANFILTERREPLICATE Running median of width k with edge replication.
%   Matches scipy.ndimage.median_filter(..., mode='nearest'), which MATLAB's
%   MEDFILT1 does not: its edge handling pads with zeros or shortens the
%   window, and either would distort the ends of a laminar profile.

x = x(:);
n = numel(x);
h = floor(k / 2);
idx = (1:n)' + (-h:h);
idx = min(max(idx, 1), n);
y = median(x(idx), 2);
end
