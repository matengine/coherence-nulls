function spectral_shape_comb_vs_bump(recordingFiles, nSeeds, durationS)
%SPECTRAL_SHAPE_COMB_VS_BUMP Can spectral shape separate a reset from a bump?
%   A generator reset on every cycle of a slow rhythm is cyclostationary at the
%   slow period, so its cross-regional coherence sits on the harmonics of that
%   rhythm. A band-limited generator gives a smooth bump. The two can share a
%   band mean, so the band mean cannot separate them. The comb index, the
%   spread of coherence across the band in units of sampling noise, is raised
%   by both; the harmonic contrast, coherence on the harmonics minus coherence
%   midway between them, is large for a comb and near zero for a bump.
%
%   Frequency drift smears a comb away, so the comparison means nothing unless
%   the slow rhythm is equally stable on both sides. Drift is measured with one
%   procedure everywhere, and the model is re-run at the drift measured in each
%   recording before the two are compared.
%
%   spectral_shape_comb_vs_bump
%   spectral_shape_comb_vs_bump({'session_csd.mat'})
%
%   A recording file is a .mat holding f, R, L and slowDriftHz, which
%   run_session_coherence writes.

if nargin < 1, recordingFiles = {}; end
if nargin < 2 || isempty(nSeeds), nSeeds = 5; end
if nargin < 3 || isempty(durationS), durationS = 300; end
addpath(fileparts(fileparts(mfilename('fullpath'))));
GAMMA = [60 100];
SWEEP = [7.5 10.0];

drifts = [0.05 0.10 0.15 0.20 0.25 0.30 0.35 0.40 0.50 0.60 0.80 1.00];
fprintf('Reset model: both statistics against the measured slow-rhythm drift\n');
fprintf('%9s %9s %16s %16s\n', 'parameter', 'measured', 'comb index', 'harmonic contrast');
rows = zeros(numel(drifts), 5);
for i = 1:numel(drifts)
    c = zeros(nSeeds, 1); h = zeros(nSeeds, 1); m = zeros(nSeeds, 1);
    for k = 1:nSeeds
        [c(k), h(k), m(k)] = runModel(drifts(i), 600 + 10 * i + k, durationS, GAMMA);
    end
    rows(i, :) = [drifts(i), mean(m), mean(c), mean(h), std(h)];
    fprintf('%9.2f %9.2f %8.2f +/-%4.2f %8.2f +/-%4.2f\n', ...
        drifts(i), mean(m), mean(c), std(c), mean(h), std(h));
end
fprintf('\nThe comb index settles on a plateau near 1.2 to 1.5 while the harmonic contrast goes to\n');
fprintf('zero. A comb index of 1.5 is equally consistent with a washed-out comb and with a modest\n');
fprintf('bump; a harmonic contrast of zero is consistent only with the bump.\n\n');

for i = 1:numel(recordingFiles)
    z = load(recordingFiles{i});
    f = z.f(:).';
    R = z.R(:).';
    L = double(z.L);
    drift = NaN;
    if isfield(z, 'slowDriftHz'), drift = z.slowDriftHz; end
    m = coherence.bandMask(f, [5 11]);
    [~, j] = max(R .* m);
    fSlow = f(j);
    [bestF, bestV] = coherence.harmonicContrastSweep(R, f, GAMMA, L, SWEEP(1), SWEEP(2));
    fprintf('%s\n', recordingFiles{i});
    fprintf('  measured drift %.2f Hz, slow peak %.1f Hz\n', drift, fSlow);
    fprintf('  comb index %.2f, harmonic contrast %+.2f (best over %.1f-%.1f Hz: %+.2f at %.1f Hz)\n', ...
        coherence.combIndex(R, f, GAMMA, L), ...
        coherence.harmonicContrast(R, f, GAMMA, L, fSlow), SWEEP(1), SWEEP(2), bestV, bestF);
    if isfinite(drift)
        par = interp1(rows(:, 2), rows(:, 1), drift, 'linear', 'extrap');
        c = zeros(8, 1); h = zeros(8, 1); mm = zeros(8, 1);
        for k = 1:8
            [c(k), h(k), mm(k)] = runModel(par, 9000 + k, durationS, GAMMA);
        end
        fprintf('  reset model at the same measured drift (%.2f achieved, 8 realisations):\n', mean(mm));
        fprintf('    comb index %.2f (%.2f-%.2f), harmonic contrast %.2f (%.2f-%.2f)\n', ...
            mean(c), min(c), max(c), mean(h), min(h), max(h));
    end
    fprintf('\n');
end
end


function [c, h, measured] = runModel(drift, seed, durationS, GAMMA)
p = models.resetParams('durationS', durationS, 'seed', seed, 'gammaJitterHz', 5, ...
    'thetaLocalAmp', 0.8, 'resetR', 0.65, 'thetaMeanHz', 8.8, ...
    'thetaDriftSdHz', drift, 'thetaDriftTauS', 2, ...
    'thetaLocalDriftSdHz', drift, 'thetaLocalDriftTauS', 2);
s = models.simulateReset(p);
measured = coherence.thetaWindowDrift(s.lfpA, p.fs);
[f, X] = coherence.segmentFFTs(s.lfpA, p.fs);
[~, Y] = coherence.segmentFFTs(s.lfpB, p.fs);
R = coherence.coherenceR(X, Y);
L = size(X, 1);
c = coherence.combIndex(R, f, GAMMA, L);
h = coherence.harmonicContrast(R, f, GAMMA, L, p.thetaMeanHz);
end
