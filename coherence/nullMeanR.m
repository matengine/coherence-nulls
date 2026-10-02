function m = nullMeanR(L)
%NULLMEANR Expected coherence under the null, (L-1)*B(3/2, L-1).
%   This is the estimator's finite-sample floor, not a property of the data.
%   The gap between an observed value and a surrogate one is therefore not an
%   additive amount of coupling.

m = (L - 1) .* beta(1.5, L - 1);
end
