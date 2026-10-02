function idx = troughsFromPhase(theta)
%TROUGHSFROMPHASE Sample indices at which a phase series passes pi upward.
%   These are the cosine minima, i.e. the troughs of the rhythm.

w = coherence.wrapAngle(theta(:) - pi);
idx = find(w(1:end - 1) < 0 & w(2:end) >= 0) + 1;
end
