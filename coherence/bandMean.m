function v = bandMean(Rf, f, band)
%BANDMEAN Mean of R over the bins of band, along the last dimension.
%   Accepts a single spectrum or a surrogates-by-frequency array.

m = coherence.bandMask(f, band);
v = mean(Rf(:, m), 2);
end
