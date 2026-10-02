function c = combIndex(Rf, f, band, L)
%COMBINDEX Spread of coherence across a band, in units of sampling noise.
%   C = COMBINDEX(RF, F, BAND, L) divides the standard deviation of R across
%   the bins of the band by the standard deviation expected from sampling alone.
%
%   One means a flat band within sampling noise. Raised by a comb and by a
%   smooth bump alike, so pair it with coherence.harmonicContrast before
%   concluding anything about the shape of a spectrum.

r = Rf(coherence.bandMask(f, band));
c = std(r, 1) / ((1 - mean(r) ^ 2) / sqrt(2 * L));
end
