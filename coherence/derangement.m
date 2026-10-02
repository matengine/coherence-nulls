function p = derangement(n)
%DERANGEMENT A permutation of 1..n with no fixed point.
%   So no segment or epoch is ever paired with itself.

if n < 2
    p = 1:n;
    return
end
while true
    p = randperm(n);
    if all(p ~= 1:n)
        return
    end
end
end
