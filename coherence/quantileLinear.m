function q = quantileLinear(x, p)
%QUANTILELINEAR Quantile along the first dimension, linear interpolation rule.
%   Q = QUANTILELINEAR(X, P) uses the same definition as numpy's default, so
%   that this port and the Python one return the same surrogate thresholds.
%   MATLAB's QUANTILE uses a different rule, which shifts a 95th percentile of
%   a few hundred surrogates by a small but visible amount.

if isvector(x)
    x = x(:);
end
xs = sort(x, 1);
n = size(xs, 1);
h = (n - 1) * p + 1;
lo = min(max(floor(h), 1), n);
hi = min(max(ceil(h), 1), n);
q = xs(lo, :) + (h - lo) .* (xs(hi, :) - xs(lo, :));
end
