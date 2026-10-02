function v = powerWeightedR2(Rf, Syy, f, band)
%POWERWEIGHTEDR2 Fraction of band power linearly predictable from the other signal.
%   Distinct from both the square of a band-mean magnitude and the mean of R^2
%   across the band. Saying which of the three is meant removes a common
%   ambiguity in reported numbers.

m = coherence.bandMask(f, band);
v = sum(Rf(:, m) .^ 2 .* Syy(:, m), 2) ./ sum(Syy(:, m), 2);
end
