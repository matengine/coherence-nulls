function w = wrapAngle(phi)
%WRAPANGLE Wrap angles in radians to (-pi, pi].

w = angle(exp(1i * phi));
end
