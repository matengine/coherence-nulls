function conditional_null_calibration(durationS, nSurr, nSeeds, couplingC)
%CONDITIONAL_NULL_CALIBRATION When does the conditional null report a false positive?
%   Three cases on the reset model, which has no coupling between the gamma
%   generators unless one is added:
%
%     c = 0, true reset boundaries    must not exceed, or the test is
%                                     anti-conservative and a positive result
%                                     means nothing
%     c > 0, true reset boundaries    must exceed, or the test has no power
%     c = 0, field-derived boundaries the realistic case for a field
%                                     recording: boundaries sit tens of
%                                     milliseconds from the reset, the
%                                     surrogate collapses, and the test fires
%                                     with no coupling present
%
%   The third case is why a positive result on field-derived boundaries cannot
%   support a claim about a channel. The failure is one-sided, so a
%   non-exceedance stays valid with any boundaries.
%
%   conditional_null_calibration
%   conditional_null_calibration(200, 200)

if nargin < 1 || isempty(durationS), durationS = 200; end
if nargin < 2 || isempty(nSurr), nSurr = 200; end
if nargin < 3 || isempty(nSeeds), nSeeds = 6; end
if nargin < 4 || isempty(couplingC), couplingC = 0.10; end
addpath(fileparts(fileparts(mfilename('fullpath'))));

labels = {'no coupling, true reset boundaries', ...
          sprintf('coupling %.2f, true reset boundaries', couplingC), ...
          'no coupling, field-derived boundaries'};
coup = [0 couplingC 0];
useField = [false false true];
expect = {'rate should be near the nominal 5%', 'should exceed', 'exceeds: FALSE POSITIVE'};

fprintf('Calibration is a rate over realisations, not a property of one run: with no coupling a\n');
fprintf('single realisation can land either side of its own 95th percentile.\n\n');
fprintf('%-38s %10s %11s %9s %14s\n', 'case', 'observed', 'surrogate95', 'exceeded', 'boundary error');
for c = 1:3
    obs = zeros(nSeeds, 1); q95 = zeros(nSeeds, 1); ex = 0; err = zeros(nSeeds, 1);
    for k = 1:nSeeds
        [obs(k), q95(k), hit, err(k)] = oneCase(coup(c), useField(c), 510 + 7 * k, durationS, nSurr);
        ex = ex + hit;
    end
    fprintf('%-38s %10.3f %11.3f %4d/%d %11.1f ms   (%s)\n', labels{c}, mean(obs), mean(q95), ...
        ex, nSeeds, mean(err), expect{c});
end
fprintf('\nRead the third line before using this test on a recording: with boundaries taken from a\n');
fprintf('field whose slow rhythm is a mixture, exceedance is what shared drive alone produces.\n');
end


function [obs, q95, hit, offsetMs] = oneCase(couplingC, useFieldBoundaries, seed, durationS, nSurr)
GAMMA = [60 100];
p = models.resetParams('durationS', durationS, 'seed', seed, 'resetR', 0.65, ...
    'thetaLocalAmp', 0.8, 'couplingC', couplingC);
s = models.simulateReset(p);
if useFieldBoundaries
    troughs = coherence.troughsFromField(s.lfpB, p.fs);
    trueT = s.troughsB;
    near = troughs(interp1(troughs, 1:numel(troughs), trueT, 'nearest', 'extrap'));
    offsetMs = median(abs(near - trueT)) / p.fs * 1000;
else
    troughs = s.troughsB;
    offsetMs = 0;
end
[obs, sur] = coherence.cyclePermutationTest(s.lfpA, s.lfpB, troughs, p.fs, GAMMA, nSurr);
q95 = coherence.quantileLinear(sur, 0.95);
hit = obs > q95;
end
