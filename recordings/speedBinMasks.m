function masks = speedBinMasks(speedSeg, ok, bins)
%SPEEDBINMASKS One logical mask per running-speed bin, excluding artifacts.
%   MASKS = SPEEDBINMASKS(SPEEDSEG, OK) uses the default bins 1-5, 5-15, 15-35
%   and above 35 cm/s. MASKS is a cell array.

if nargin < 3 || isempty(bins)
    bins = [1 5; 5 15; 15 35; 35 Inf];
end
masks = cell(size(bins, 1), 1);
for i = 1:size(bins, 1)
    masks{i} = ok(:) & speedSeg(:) > bins(i, 1) & speedSeg(:) <= bins(i, 2);
end
end
