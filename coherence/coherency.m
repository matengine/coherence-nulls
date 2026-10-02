function C = coherency(X, Y)
%COHERENCY Complex coherency between two spectral arrays.
%   C = COHERENCY(X, Y) averages the cross-spectrum over segments and
%   normalises by both auto-spectra.
%
%   Spectra are pooled before normalising. Normalising each segment first and
%   averaging afterwards gives a different and generally larger quantity, which
%   is one of the ways two papers reporting "coherence" disagree numerically
%   while using the same word for it.

Sxy = mean(X .* conj(Y), 1);
Sxx = mean(abs(X) .^ 2, 1);
Syy = mean(abs(Y) .^ 2, 1);
C = Sxy ./ sqrt(Sxx .* Syy);
end
