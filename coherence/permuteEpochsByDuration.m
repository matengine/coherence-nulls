function out = permuteEpochsByDuration(E, durs)
%PERMUTEEPOCHSBYDURATION Derange cycle identity within groups of equal duration.
%   Every slot receives an epoch of exactly its own length and taper, so the
%   permutation changes which cycle is paired with which and nothing else.
%   Epochs with no duration partner stay where they are.

out = E;
u = unique(durs);
for i = 1:numel(u)
    idx = find(durs == u(i));
    if numel(idx) >= 2
        out(idx, :) = E(idx(coherence.derangement(numel(idx))), :);
    end
end
end
