function prof = laminarProfile(X, fs, speed, opts)
%LAMINARPROFILE Per-channel quantities from which sites are assigned to layers.
%   PROF = LAMINARPROFILE(X, FS) with X as (channel, sample) in depth order.
%   PROF = LAMINARPROFILE(X, FS, SPEED) separates running from rest using a
%   per-sample speed. PROF = LAMINARPROFILE(X, FS, SPEED, OPTS) overrides the
%   defaults in the struct OPTS.
%
%   Returns band power during running and rest, the theta phase profile whose
%   reversal marks the boundary between two current generators, the
%   theta-trough-triggered average and its current source density, and the
%   ripple-triggered average.
%
%   These are physiological inferences, not anatomy. A single live channel can
%   inflate every laminar contrast metric, so read the whole profile rather
%   than a summary number, and check the returned correlation matrix to
%   confirm that channel order follows depth.
%
%   Fields: corr, f, Prun, Prest, thetaPower, gammaPower, ripplePower,
%   thetaPhase, cohWithRef, thetaTroughAvg, thetaTroughCsd, rippleAvg, ref,
%   rippleCh, thetaPeakHz, nTroughs, nRipples, nRunSamples, nRestSamples,
%   restIsFallback.

if nargin < 3, speed = []; end
d = struct('runAbove', 15, 'restBelow', 5, 'thetaBand', [6 10], 'rippleBand', [150 250], ...
           'gammaBand', [60 100], 'troughHalfwidthS', 0.075, 'rippleHalfwidthS', 0.05, 'rippleSd', 5);
if nargin >= 4
    fn = fieldnames(opts);
    for i = 1:numel(fn), d.(fn{i}) = opts.(fn{i}); end
end

X = double(X);
X = X - mean(X, 2);
[nch, nsamp] = size(X);
minSamples = round(2 * fs);
restIsFallback = false;
if isempty(speed)
    run = true(1, nsamp);
    rest = true(1, nsamp);
else
    v = speed(:).';
    finite = isfinite(v);
    run = finite & v > d.runAbove;
    rest = finite & v < d.restBelow;
    % A session can contain almost no running, or almost no rest, and a
    % spectral estimate on an empty selection is an error rather than a
    % warning. Fall back to the fastest and slowest fifths, then to everything.
    if sum(run) < minSamples && any(finite)
        run = finite & v >= quantile(v(finite), 0.8);
    end
    if sum(run) < minSamples, run = true(1, nsamp); end
    restIsFallback = sum(rest) < minSamples;
    if restIsFallback && any(finite)
        rest = finite & v <= quantile(v(finite), 0.2);
    end
    if sum(rest) < minSamples, rest = true(1, nsamp); end
end

[b, a] = butter(2, [1 100] / (fs / 2), 'bandpass');
C = corrcoef(filtfilt(b, a, X.'));

nper = round(fs);
[Prun, f] = pwelch(X(:, run).', coherence.hannPeriodic(nper), floor(nper / 2), nper, fs);
Prest = pwelch(X(:, rest).', coherence.hannPeriodic(nper), floor(nper / 2), nper, fs);
Prun = Prun.'; Prest = Prest.'; f = f.';
th = f >= d.thetaBand(1) & f <= d.thetaBand(2);
rip = f >= d.rippleBand(1) & f <= d.rippleBand(2);
gam = f >= d.gammaBand(1) & f <= d.gammaBand(2);
thetaPower = mean(Prun(:, th), 2);
gammaPower = mean(Prun(:, gam), 2);
ripplePower = mean(Prest(:, rip), 2);

[~, ref] = max(thetaPower);
[Pxy, fc] = cpsd(repmat(X(ref, run).', 1, nch), X(:, run).', coherence.hannPeriodic(nper), floor(nper / 2), nper, fs);
[~, iTh] = max(Prun(ref, th));
fTh = f(find(th, 1) + iTh - 1);
[~, ipk] = min(abs(fc - fTh));
thetaPhase = rad2deg(angle(Pxy(ipk, :))).';
cohWithRef = (abs(Pxy(ipk, :)).' ./ sqrt(Prun(ref, ipk) * Prun(:, ipk)));

[bt, at] = butter(2, [5 11] / (fs / 2), 'bandpass');
ph = angle(hilbert(filtfilt(bt, at, X(ref, :).')));
w = coherence.wrapAngle(ph - pi);
troughs = find(w(1:end - 1) < 0 & w(2:end) >= 0) + 1;
troughs = troughs(run(troughs).');
half = round(d.troughHalfwidthS * fs);
troughs = troughs(troughs > half & troughs < nsamp - half);
ta = zeros(nch, 2 * half);
for k = 1:numel(troughs)
    ta = ta + X(:, troughs(k) - half:troughs(k) + half - 1);
end
ta = ta / max(numel(troughs), 1);
csd = -diff(ta, 2, 1);

[~, rch] = max(ripplePower);
[br, ar] = butter(2, d.rippleBand / (fs / 2), 'bandpass');
env = abs(hilbert(filtfilt(br, ar, X(rch, :).')));
thr = mean(env(rest)) + d.rippleSd * std(env(rest));
gated = env .* rest.';
if any(gated > thr)
    [~, peaks] = findpeaks(gated, 'MinPeakHeight', thr, 'MinPeakDistance', round(0.05 * fs));
else
    peaks = [];        % no ripple-like events, which is normal in some recordings
end
halfr = round(d.rippleHalfwidthS * fs);
peaks = peaks(peaks > halfr & peaks < nsamp - halfr);
ra = zeros(nch, 2 * halfr);
for k = 1:numel(peaks)
    ra = ra + X(:, peaks(k) - halfr:peaks(k) + halfr - 1);
end
ra = ra / max(numel(peaks), 1);

prof = struct('corr', C, 'f', f, 'Prun', Prun, 'Prest', Prest, ...
    'thetaPower', thetaPower, 'gammaPower', gammaPower, 'ripplePower', ripplePower, ...
    'thetaPhase', thetaPhase, 'cohWithRef', cohWithRef, ...
    'thetaTroughAvg', ta, 'thetaTroughCsd', csd, 'rippleAvg', ra, ...
    'ref', ref, 'rippleCh', rch, 'thetaPeakHz', fTh, ...
    'nTroughs', numel(troughs), 'nRipples', numel(peaks), ...
    'nRunSamples', sum(run), 'nRestSamples', sum(rest), 'restIsFallback', restIsFallback);
end
