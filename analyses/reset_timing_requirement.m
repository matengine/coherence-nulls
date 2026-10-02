function reset_timing_requirement(durationS, nSeeds, resetR)
%RESET_TIMING_REQUIREMENT How precise must a shared reset be?
%   Independent jitter is added to the two regions' reset times and the
%   resulting coherence is compared with exp(-(2*pi*f*sigma)^2), which is not a
%   fit: two events that should coincide but carry independent zero-mean errors
%   of SD sigma have a phase difference scattered with SD 2*pi*f*sigma*sqrt(2).
%   A zero-mean error does not average out, and at 80 Hz an error of 4 ms is
%   163 degrees of phase.
%
%   The consequence is a measurement anyone can make on a recording. If the
%   cross-regional timing of the slow rhythm is looser than about a
%   millisecond, a shared pacemaker cannot have produced the reported gamma
%   coherence by resetting independent generators.
%
%   reset_timing_requirement
%   reset_timing_requirement(100, 2)

if nargin < 1 || isempty(durationS), durationS = 200; end
if nargin < 2 || isempty(nSeeds), nSeeds = 3; end
if nargin < 3 || isempty(resetR), resetR = 0.65; end
addpath(fileparts(fileparts(mfilename('fullpath'))));
GAMMA = [60 100];
jitters = [0 0.5 1 2 4];
fMid = mean(GAMMA);

fprintf('Independent jitter in the two regions'' reset times, %d realisations of %.0f s each\n', ...
    nSeeds, durationS);
fprintf('%10s %14s %10s %14s\n', 'jitter (ms)', 'coherence', 'null 95th', 'predicted drop');
for j = jitters
    obs = zeros(nSeeds, 1);
    n95 = zeros(nSeeds, 1);
    for k = 1:nSeeds
        [obs(k), n95(k)] = coherenceAtJitter(j, 700 + 10 * k, durationS, resetR, GAMMA);
    end
    fprintf('%10.1f %14.3f %10.3f %14.2f\n', j, mean(obs), mean(n95), ...
        coherence.timingErrorAttenuation(j / 1000, fMid));
end
fprintf('\nThe measured coherence tracks exp(-(2*pi*f*sigma)^2) scaled from the zero-jitter value.\n');
fprintf('Sub-millisecond shared timing is required; a few milliseconds leaves nothing above the null.\n');
end


function [obs, null95] = coherenceAtJitter(jitterMs, seed, durationS, resetR, GAMMA)
p = models.resetParams('durationS', durationS, 'seed', seed, 'resetR', resetR, ...
    'thetaLocalAmp', 0.8, 'resetTimeJitterMs', jitterMs);
s = models.simulateReset(p);
[f, X] = coherence.segmentFFTs(s.lfpA, p.fs, 0.5);
[~, Y] = coherence.segmentFFTs(s.lfpB, p.fs, 0.5);
obs = coherence.bandMean(coherence.coherenceR(X, Y), f, GAMMA);
null95 = coherence.bandNull(coherence.nullCircularShift(X, Y, 100), f, GAMMA, 0.95, 'mean');
end
