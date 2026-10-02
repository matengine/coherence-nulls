function [f, X] = rfftRows(xs, fs)
%RFFTROWS One-sided FFT of each row, with its frequency axis.
%   [F, X] = RFFTROWS(XS, FS) returns the non-negative-frequency half of the
%   row-wise FFT and the matching frequency vector.

n = size(xs, 2);
F = fft(xs, [], 2);
X = F(:, 1:(floor(n / 2) + 1));
f = (0:floor(n / 2)) * fs / n;
end
