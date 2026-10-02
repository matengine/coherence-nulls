function P = autoSpectrum(X)
%AUTOSPECTRUM Mean power spectrum over the segments of a spectral array.

P = mean(abs(X) .^ 2, 1);
end
