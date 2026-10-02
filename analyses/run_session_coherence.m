function run_session_coherence(blockDir, probeA, siteA, probeB, siteB, nChannels, nSurr, fastBin, outPrefix, maxSeconds)
%RUN_SESSION_COHERENCE Cross-regional coherence of one session, by speed bin.
%   RUN_SESSION_COHERENCE(BLOCKDIR, PROBEA, SITEA, PROBEB, SITEB) runs the
%   whole pipeline: gain equalisation, bad-channel detection, the current
%   source density, segmentation by speed, coherence on both the field
%   potential and the current source density, the circular-shift and
%   segment-permutation nulls, the frequency-specific effective sample size,
%   the slow-rhythm drift, and the cross-regional trough timing.
%
%   Choose the two sites with run_laminar_profile first.
%
%   Writes one .mat per derivation holding f, R, L and slowDriftHz for the
%   fastest bin, which is what spectral_shape_comb_vs_bump reads.
%
%   run_session_coherence('F:\tank\779-180601-112347', 'Wav4', 20, 'Wav3', 16)
%
%   Pass MAXSECONDS to read only the first part of a session, which is useful
%   for a quick look before committing to the full read.

if nargin < 6 || isempty(nChannels), nChannels = 32; end
if nargin < 7 || isempty(nSurr), nSurr = 200; end
if nargin < 8 || isempty(fastBin), fastBin = 4; end
if nargin < 9 || isempty(outPrefix), outPrefix = 'session'; end
if nargin < 10, maxSeconds = []; end
addpath(fileparts(fileparts(mfilename('fullpath'))));
rng(7, 'twister');
THETA = [5 11];
binLabels = {'1-5', '5-15', '15-35', '>35'};

d = recordings.readStreams(blockDir, {probeA, probeB}, 1:nChannels, 0, maxSeconds, [], [], [], true);
fs = d.fs;
A = d.data.(probeA);
B = d.data.(probeB);
nsamp = min(size(A, 2), size(B, 2));
A = A(:, 1:nsamp); B = B(:, 1:nsamp);

tr = recordings.readTracking(blockDir);
[tg, v] = recordings.speedFromTracking(tr);
speed = interp1(tg, v, (0:nsamp - 1) / fs, 'linear', NaN);

Aeq = recordings.equalizeGains(A, fs);
Beq = recordings.equalizeGains(B, fs);
badA = recordings.badChannels(Aeq, fs);
badB = recordings.badChannels(Beq, fs);
fprintf('bad channels: %s %s, %s %s\n', probeA, mat2str(badA), probeB, mat2str(badB));
Acsd = recordings.currentSourceDensity(Aeq, badA);
Bcsd = recordings.currentSourceDensity(Beq, badB);

kinds = {'lfp', 'csd'};
traces = {{A(siteA, :), B(siteB, :)}, {Acsd(siteA, :), Bcsd(siteB, :)}};
for ki = 1:2
    xa = traces{ki}{1};
    xb = traces{ki}{2};
    [f, F, speedSeg, ok] = recordings.segmentSession(struct('a', xa, 'b', xb), fs, speed);
    masks = recordings.speedBinMasks(speedSeg, ok);
    fprintf('\n%s, segments per speed bin: %s\n', upper(kinds{ki}), ...
        mat2str(cellfun(@sum, masks).'));
    fprintf('%7s %6s %9s %20s %26s\n', 'bin', 'L', 'L_eff(th)', 'theta peak (null95)', ...
        'gamma mean (null95), sigma');
    store = [];
    for bi = 1:numel(masks)
        if sum(masks{bi}) < 20, continue; end
        r = recordings.pairCoherence(F.a, F.b, masks{bi}, f, nSurr);
        th = r.bands.theta;
        ga = r.bands.gamma;
        fprintf('%7s %6d %9.0f %10.3f (%.3f)      %8.3f (%.3f)  %5.0f deg\n', ...
            binLabels{bi}, r.L, th.Leff, th.observed, th.null95, ga.observed, ga.null95, ga.sigmaDeg);
        if bi == fastBin
            [drift, meanHz, nwin] = coherence.thetaWindowDrift(xb, fs, THETA, 1.0, masks{bi});
            store = struct('f', f, 'R', r.R, 'L', r.L, 'Leff', r.Leff, ...
                'q95shift', r.q95shift, 'q95perm', r.q95perm, ...
                'slowDriftHz', drift, 'slowMeanHz', meanHz, 'slowNWindows', nwin, ...
                'thetaObserved', th.observed, 'thetaNull95', th.null95, ...
                'gammaObserved', ga.observed, 'gammaNull95', ga.null95);
            fprintf('        slow-rhythm drift in this bin: %.2f Hz about %.2f Hz over %d windows\n', ...
                drift, meanHz, nwin);
        end
    end
    if ~isempty(store)
        path = sprintf('%s_%s.mat', outPrefix, kinds{ki});
        save(path, '-struct', 'store');
        fprintf('        wrote %s\n', path);
    end

    fast = false(nsamp, 1);
    n = round(fs);
    for k = find(masks{fastBin}).'
        fast((k - 1) * n + 1:min(k * n, nsamp)) = true;
    end
    off = recordings.troughOffsets(xa, xb, fs, THETA, fast);
    fprintf('  trough offsets, fastest bin: n %d, mean %.1f ms, SD %.1f ms, circular SD %.0f deg\n', ...
        off.n, off.meanMs, off.sdMs, off.circularSdDeg);
end

fprintf('\nIf the field-potential coherence has no spectral peak and climbs toward the Nyquist band\n');
fprintf('while the current source density at the same sites does not, the field-potential value is\n');
fprintf('a shared broadband component and is not the quantity to interpret.\n');
end
