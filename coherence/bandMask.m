function m = bandMask(f, band)
%BANDMASK Logical mask of the frequency bins inside band, endpoints included.

m = f >= band(1) & f <= band(2);
end
