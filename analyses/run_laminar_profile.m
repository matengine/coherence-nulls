function run_laminar_profile(blockDir, probes, nChannels, windowS, outFile)
%RUN_LAMINAR_PROFILE Laminar profile of both probes of one session.
%   RUN_LAMINAR_PROFILE(BLOCKDIR) reads a TDT tank, takes the window with the
%   most running, and saves the per-channel quantities from which layers are
%   conventionally identified, together with candidate sites.
%
%   Nothing is plotted. Inspect the saved arrays, and look at the whole profile
%   rather than at a single summary number: one live channel in a dead probe
%   can inflate every laminar contrast metric.
%
%   For another acquisition system, replace the reading step and call
%   recordings.laminarProfile on your own (channel, sample) array.
%
%   run_laminar_profile('F:\tank\779-180601-112347')

if nargin < 2 || isempty(probes), probes = {'Wav3', 'Wav4'}; end
if nargin < 3 || isempty(nChannels), nChannels = 32; end
if nargin < 4 || isempty(windowS), windowS = 600; end
if nargin < 5 || isempty(outFile), outFile = 'laminar_profile.mat'; end
addpath(fileparts(fileparts(mfilename('fullpath'))));

tr = recordings.readTracking(blockDir);
[tg, v] = recordings.speedFromTracking(tr);
[t0, frac] = bestWindow(tg, v, windowS);
fprintf('window %.0f to %.0f s, %.0f%% of it above 15 cm/s\n', t0, t0 + windowS, 100 * frac);

d = recordings.readStreams(blockDir, probes, 1:nChannels, t0, t0 + windowS, [], [], [], false);
fs = d.fs;
saved = struct('windowStartS', t0, 'windowS', windowS, 'fs', fs);
for i = 1:numel(probes)
    store = probes{i};
    X = d.data.(store);
    speed = interp1(tg, v, t0 + (0:size(X, 2) - 1) / fs, 'linear', 'extrap');
    prof = recordings.laminarProfile(X, fs, speed);
    sites = recordings.suggestSites(prof);
    fprintf('\n%s: theta reference channel %d at %.1f Hz, ripple-power maximum channel %d, ', ...
        store, prof.ref, prof.thetaPeakHz, prof.rippleCh);
    fprintf('%d troughs, %d ripple events\n', prof.nTroughs, prof.nRipples);
    fprintf('  median correlation with the next channel %.2f, eight away %.2f  ', ...
        median(diag(prof.corr, 1)), median(diag(prof.corr, 8)));
    fprintf('(a low second value means channel order follows depth)\n');
    span = max(prof.thetaPhase) - min(prof.thetaPhase);
    if isnan(sites.phaseReversal)
        fprintf('  theta phase spans %.0f degrees; no reversal, so no phase landmark\n', span);
    else
        fprintf('  theta phase spans %.0f degrees, reversal near channel %d\n', span, sites.phaseReversal);
    end
    fprintf('  candidate sites: pyramidal %d, theta maximum %d, theta sink %d\n', ...
        sites.pyramidal, sites.thetaMax, sites.thetaSink);
    saved.(store) = prof;
    saved.([store '_sites']) = sites;
end

save(outFile, '-struct', 'saved');
fprintf('\nwrote %s\n', outFile);
fprintf('These are physiological inferences. Call a site a layer only with histology registered\n');
fprintf('to the probe geometry.\n');
end


function [t0, frac] = bestWindow(t, v, winS, above, stepS)
if nargin < 4, above = 15; end
if nargin < 5, stepS = 30; end
best = -1; t0 = 0;
for ts = 0:stepS:max(t(end) - winS, 1)
    m = t >= ts & t < ts + winS;
    if ~any(m), continue; end
    fr = mean(v(m) > above);
    if fr > best, best = fr; t0 = ts; end
end
frac = best;
end
