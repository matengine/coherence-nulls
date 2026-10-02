function L = lForSignificance(R, alpha)
%LFORSIGNIFICANCE Smallest L at which an observed R clears the pointwise test.

if nargin < 2, alpha = 0.05; end
L = 1 + log(alpha) ./ log(1 - R .^ 2);
end
