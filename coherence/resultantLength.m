function R = resultantLength(phases)
%RESULTANTLENGTH Mean resultant length of phases in radians, along dimension 1.

R = abs(mean(exp(1i * phases), 1));
end
