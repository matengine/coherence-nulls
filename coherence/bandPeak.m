function v = bandPeak(Rf, f, band)
%BANDPEAK Maximum of R over the bins of band, along the last dimension.
%   Noisier than the mean, with a correspondingly higher null, so it needs its
%   own surrogate rather than the mean's.

m = coherence.bandMask(f, band);
v = max(Rf(:, m), [], 2);
end
