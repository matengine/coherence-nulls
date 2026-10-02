function detection_null_calibration(nPairs, durationS, segmentS, nSurr, atHz)
%DETECTION_NULL_CALIBRATION Does the detection null have the rate it claims?
%   Why a detection null has to preserve each signal's serial structure, and
%   what persistence costs at the frequency of a rhythm.
%
%   On independent model recordings, where the true coupling is zero, this
%   compares five nulls at a few frequencies: the empirical one from many
%   independent pairs, which is the truth; the segment-permutation null, which
%   sits below it at the rhythm's frequency because re-pairing destroys the
%   persistence the real data have; the circular-shift null, which tracks it;
%   and the analytic null evaluated at the raw segment count and at L_eff(f).
%   It then reports the false-positive rate each surrogate delivers at a
%   nominal 5%.
%
%   detection_null_calibration
%   detection_null_calibration(40, 100)

if nargin < 1 || isempty(nPairs), nPairs = 30; end
if nargin < 2 || isempty(durationS), durationS = 100; end
if nargin < 3 || isempty(segmentS), segmentS = 0.5; end
if nargin < 4 || isempty(nSurr), nSurr = 200; end
if nargin < 5 || isempty(atHz), atHz = [8 80]; end
addpath(fileparts(fileparts(mfilename('fullpath'))));
rng(3, 'twister');

trueR = [];
for k = 1:nPairs
    [f, X, Y] = independentPair(2000 + k, durationS, segmentS);
    trueR(k, :) = coherence.coherenceR(X, Y); %#ok<AGROW>
end

[f, X, Y] = independentPair(2000, durationS, segmentS);
N = size(X, 1);
Rshift = coherence.nullCircularShift(X, Y, nSurr);
Rperm = coherence.nullPermutation(X, Y, nSurr);
Leff = coherence.lEffFrequency(X, Y, min(30, floor(N / 4)));

fprintf('Independent recordings, %d segments of %.1f s, %d pairs for the true null\n\n', ...
    N, segmentS, nPairs);
fprintf('%8s %10s %10s %10s %12s %12s %10s\n', ...
    'freq', 'true 95th', 'shift 95th', 'perm 95th', 'analytic L=N', 'analytic Leff', 'Leff');
for t = atHz
    [~, i] = min(abs(f - t));
    fprintf('%7.1f %10.3f %10.3f %10.3f %12.3f %12.3f %10.0f\n', f(i), ...
        coherence.quantileLinear(trueR(:, i), 0.95), ...
        coherence.quantileLinear(Rshift(:, i), 0.95), ...
        coherence.quantileLinear(Rperm(:, i), 0.95), ...
        coherence.critR(N), coherence.critR(max(Leff(i), 2)), Leff(i));
end

fprintf('\nFalse-positive rate at a nominal 5%%, over the %d independent pairs:\n', nPairs);
for t = atHz
    [~, i] = min(abs(f - t));
    fpPerm = mean(trueR(:, i) > coherence.quantileLinear(Rperm(:, i), 0.95));
    fpShift = mean(trueR(:, i) > coherence.quantileLinear(Rshift(:, i), 0.95));
    fpAnal = mean(trueR(:, i) > coherence.critR(N));
    fprintf('  %5.1f Hz   permutation %4.0f%%   circular shift %4.0f%%   analytic at L=N %4.0f%%\n', ...
        f(i), 100 * fpPerm, 100 * fpShift, 100 * fpAnal);
end
fprintf('\nPersistence is frequency-specific: the effective count collapses at the rhythm''s own\n');
fprintf('frequency and is close to the segment count elsewhere, so a flat threshold is wrong\n');
fprintf('exactly in the band where a rhythm lives.\n');
end


function [f, X, Y] = independentPair(seed, durationS, segmentS)
% Two recordings that share nothing: separate simulations with different seeds.
a = models.simulateReset(models.resetParams('durationS', durationS, 'seed', seed, 'thetaLocalAmp', 0.8));
b = models.simulateReset(models.resetParams('durationS', durationS, 'seed', seed + 5000, 'thetaLocalAmp', 0.8));
[f, X] = coherence.segmentFFTs(a.lfpA, a.fs, segmentS);
[~, Y] = coherence.segmentFFTs(b.lfpB, b.fs, segmentS);
end
