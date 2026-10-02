function Y = bipolar(X)
%BIPOLAR Adjacent-channel difference, returning n-1 traces.
%   The other derivation that removes a shared broadband component. Use it
%   where the probe geometry does not justify a second difference, for example
%   on tetrodes or where channel spacing is not uniform.

Y = diff(double(X), 1, 1);
end
