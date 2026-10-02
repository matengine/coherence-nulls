function sites = suggestSites(prof)
%SUGGESTSITES Candidate recording sites from a laminar profile. Inferences only.
%   SITES = SUGGESTSITES(PROF) returns a struct of channel indices:
%
%     pyramidal     maximum ripple-band power during rest
%     thetaMax      maximum theta-band power during running
%     thetaSink     largest sink in the theta-trough-triggered CSD, the local
%                   dendritic input zone. Use this in entorhinal cortex, where
%                   the theta phase profile is often flat and gives no landmark
%     phaseReversal where the theta phase profile crosses the midpoint of its
%                   range, or NaN if the profile spans under 90 degrees
%
%   Calling any of these a layer requires histology registered to the probe
%   geometry. Without that, say "theta-sink site", not "layer III".

[~, pyr] = max(prof.ripplePower);
[~, thMax] = max(prof.thetaPower);
[~, linIdx] = min(prof.thetaTroughCsd(:));
[row, ~] = ind2sub(size(prof.thetaTroughCsd), linIdx);
sink = row + 1;                       % csd row i corresponds to channel i+1
phase = prof.thetaPhase;
rev = NaN;
if max(phase) - min(phase) > 90
    target = 0.5 * (max(phase) + min(phase));
    [~, rev] = min(abs(phase - target));
end
sites = struct('pyramidal', pyr, 'thetaMax', thMax, 'thetaSink', sink, 'phaseReversal', rev);
end
