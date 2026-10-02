function s = sigmaFromR(R)
%SIGMAFROMR Constant-amplitude-equivalent circular SD in degrees, sqrt(-2 ln R).
%   This is what makes a magnitude interpretable: R = 0.25 is about 95 degrees
%   of phase scatter, which constrains any single cycle only loosely however
%   reliable the population mean may be.

s = rad2deg(sqrt(-2 * log(R)));
end
