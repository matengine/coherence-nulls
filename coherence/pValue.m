function p = pValue(R, L)
%PVALUE Pointwise p-value under independence, (1 - R^2)^(L-1).

p = (1 - R .^ 2) .^ (L - 1);
end
