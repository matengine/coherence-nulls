function a = timingErrorAttenuation(sigmaS, fHz)
%TIMINGERRORATTENUATION How far independent timing error knocks down coherence.
%   A = TIMINGERRORATTENUATION(SIGMAS, FHZ) returns exp(-(2*pi*f*sigma)^2).
%
%   Two events that should coincide but carry independent zero-mean errors of
%   SD sigma have a phase difference scattered with SD 2*pi*f*sigma*sqrt(2), so
%   the resultant length falls by this factor. A zero-mean error does not
%   average out: what matters is its spread against the period. At 80 Hz the
%   factor is 0.78 at 1 ms, 0.36 at 2 ms and 0.02 at 4 ms.

a = exp(-(2 * pi * fHz .* sigmaS) .^ 2);
end
