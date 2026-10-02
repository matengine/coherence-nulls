% RECORDINGS  Reading a laminar probe recording and analysing it.
%
% Reading. Only this part is acquisition-specific; everything after it takes
% plain (channel, sample) arrays, so another system needs only a new reader.
%   readStreams            - TDT tank stream stores, resampled, without the SDK
%   readTracking           - the two-LED tracking log
%   speedFromTracking      - running speed in cm/s from the LED midpoint
%
% Preparing the signal. The spatial derivation is the step that decides what a
% cross-regional coherence means: a field-potential coherence with no spectral
% peak that climbs toward the Nyquist band is a shared broadband component, not
% a channel between regions.
%   equalizeGains          - per-channel gain from the high-frequency noise floor
%   badChannels            - dead, noisy, or neighbour-unrelated channels
%   currentSourceDensity   - negative second spatial difference
%   bipolar                - adjacent-channel difference
%   commonAverageReference - subtract the probe mean
%
% Choosing sites. These are physiological inferences, not anatomy.
%   laminarProfile         - the per-channel quantities layers are read from
%   suggestSites           - candidate sites read off a profile
%
% The cross-regional analysis
%   segmentSession         - segment by running speed, with an artifact mask
%   speedBinMasks          - one mask per speed bin
%   pairCoherence          - coherence of one pair in one bin, with its nulls
%   troughOffsets          - cross-regional timing of the slow rhythm
%
% Internal
%   analysisSlice          - a valid time window, falling back to the whole run
%   medianFilterReplicate  - running median with edge replication
