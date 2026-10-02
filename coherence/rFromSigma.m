function R = rFromSigma(sigmaDeg)
%RFROMSIGMA Resultant length implied by a circular SD in degrees.

s = deg2rad(sigmaDeg);
R = exp(-s .^ 2 / 2);
end
