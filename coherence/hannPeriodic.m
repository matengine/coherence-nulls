function w = hannPeriodic(n)
%HANNPERIODIC Periodic Hann window of length n, as a column vector.
%   A Fourier analysis wants the periodic window, which is what scipy returns
%   from get_window(..., fftbins=true). MATLAB's HANN defaults to the symmetric
%   window, which is a different one, so this is spelled out rather than left to
%   a toolbox default.

w = 0.5 * (1 - cos(2 * pi * (0:n - 1)' / n));
end
