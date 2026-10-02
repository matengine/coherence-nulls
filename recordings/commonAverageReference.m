function Y = commonAverageReference(X)
%COMMONAVERAGEREFERENCE Subtract the across-channel mean of the same probe.
%   The simplest test of whether a high-frequency coherence is common-mode: if
%   it survives this, it is not.

X = double(X);
Y = X - mean(X, 1);
end
