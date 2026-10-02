% COHERENCE  Estimating a cross-regional coherence and saying what it establishes.
%
% Spectra and the estimator
%   segmentFFTs             - non-overlapping Hann-tapered segments, as a spectral array
%   segmentFFTsMultitaper   - Slepian multitaper version
%   coherency               - complex coherency, spectra pooled before normalising
%   coherenceR              - unsquared magnitude coherence R
%   autoSpectrum            - mean power spectrum over segments
%
% How many independent estimates entered it
%   laggedSelfCoherence     - how much of a segment's spectrum the next ones repeat
%   lEffFrequency           - effective estimate count per frequency
%   lEffWithinBin           - the same when only a subset of segments is analysed
%
% What a value means at that count
%   critR                   - pointwise critical coherence
%   pValue                  - pointwise p-value under independence
%   nullMeanR               - the estimator's finite-sample floor
%   lForSignificance        - smallest L at which an observed R clears the test
%   biasCorrectedR2         - approximately unbiased population R squared
%
% What it means as phase
%   sigmaFromR              - circular SD in degrees from a resultant length
%   rFromSigma              - the inverse
%   resultantLength         - mean resultant length of a set of phases
%   wrappedNormalPdf        - density of a wrapped normal
%   timingErrorAttenuation  - how far independent timing error knocks coherence down
%   wrapAngle               - wrap radians to (-pi, pi]
%   troughsFromPhase        - indices of the cosine minima of a phase series
%
% Is it more than chance
%   derangement             - permutation with no fixed point
%   nullPermutation         - re-pairing null, anti-conservative, kept for comparison
%   nullCircularShift       - circular-shift null, the one to use
%   nullCircularShiftSubset - the same when only a subset is analysed
%   bandNull                - null distribution of a band statistic
%
% Band summaries and the shape of a spectrum
%   bandMask                - bins inside a band
%   bandMean                - band average
%   bandPeak                - band maximum
%   powerWeightedR2         - predictable fraction of band power
%   combIndex               - spread across the band in units of sampling noise
%   harmonicContrast        - coherence on the harmonics minus between them
%   harmonicContrastSweep   - its largest value over a sweep of the fundamental
%   thetaWindowDrift        - how stable the slow rhythm's frequency is
%
% Is it more than shared rhythmic drive
%   troughsFromField        - cycle boundaries from a band-passed field
%   cycleEpochFFTs          - one spectral estimate per cycle
%   permuteEpochsByDuration - derange cycle identity within equal durations
%   cyclePermutationTest    - the conditional null, with its one-sided limit
%
% Keeping this port aligned with the Python one
%   quantileLinear          - numpy's quantile rule
%   hannPeriodic            - the periodic Hann window
%   rfftRows                - one-sided row-wise FFT
