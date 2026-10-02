function out = readStreams(blockDir, stores, channels, tStart, tStop, fsOut, subBlockS, overlapS, verbose)
%READSTREAMS Read TDT tank stream stores without the TDT SDK, and resample.
%   OUT = READSTREAMS(BLOCKDIR, STORES, CHANNELS, TSTART, TSTOP) returns a
%   struct with fs, t0, fsRaw and data, where data has one field per store
%   holding a (channel, sample) array.
%
%   STORES is a cell array of four-character store names such as {'Wav3'}.
%   Pass TSTOP empty to read to the end.
%
%   The .tev payload is read in contiguous byte ranges covering all requested
%   channels at once, so an external disk is read sequentially rather than
%   seek by seek. That is the difference between minutes and hours on a
%   spinning drive.
%
%   This is the only acquisition-specific part of the package. For another
%   system, write your own reader; everything downstream takes plain
%   (channel, sample) arrays.

if nargin < 4 || isempty(tStart), tStart = 0; end
if nargin < 5, tStop = []; end
if nargin < 6 || isempty(fsOut), fsOut = 1250; end
if nargin < 7 || isempty(subBlockS), subBlockS = 60; end
if nargin < 8 || isempty(overlapS), overlapS = 1; end
if nargin < 9 || isempty(verbose), verbose = true; end

[tsqPath, tevPath] = findFiles(blockDir);
rec = readIndex(tsqPath);
startMask = rec.type == hex2dec('8801');
if any(startMask)
    t0 = rec.timestamp(find(startMask, 1));
else
    pos = rec.timestamp(rec.timestamp > 0);
    t0 = min(pos);
end

nStore = numel(stores);
idx = cell(1, nStore);
rel = cell(1, nStore);
for s = 1:nStore
    code = typecast(uint8(stores{s}), 'uint32');
    m = rec.code == code & rec.type == hex2dec('8101');
    idx{s} = structfun(@(v) v(m), rec, 'UniformOutput', false);
    rel{s} = idx{s}.timestamp - t0;
end

fsRaw = median(idx{1}.frequency);
nper = median(idx{1}.size) - 10;
fmt = formatName(median(idx{1}.format));
itemBytes = formatBytes(median(idx{1}.format));
chunkBytes = nper * itemBytes;
if isempty(tStop)
    tStop = max(cellfun(@max, rel)) + nper / fsRaw;
end
up = 64;
down = round(fsRaw * 64 / fsOut);

pieces = cell(nStore, numel(channels));
fid = fopen(tevPath, 'r');
cleaner = onCleanup(@() fclose(fid));
for b0 = tStart:subBlockS:tStop
    b1 = min(b0 + subBlockS, tStop);
    lo = b0 - overlapS;
    hi = b1 + overlapS;
    sel = cell(1, nStore);
    offs = [];
    for s = 1:nStore
        t = rel{s};
        m = (t + nper / fsRaw > lo) & (t < hi) & ismember(idx{s}.channel, channels);
        sel{s} = find(m);
        offs = [offs; idx{s}.offset(sel{s})]; %#ok<AGROW>
    end
    if isempty(offs)
        continue
    end
    oMin = min(offs);
    oMax = max(offs) + chunkBytes;
    fseek(fid, oMin, 'bof');
    buf = fread(fid, oMax - oMin, '*uint8');
    for s = 1:nStore
        t = rel{s};
        for ci = 1:numel(channels)
            ii = sel{s}(idx{s}.channel(sel{s}) == channels(ci));
            if isempty(ii)
                continue
            end
            [~, order] = sort(t(ii));
            ii = ii(order);
            x = zeros(nper, numel(ii));
            for k = 1:numel(ii)
                a = idx{s}.offset(ii(k)) - oMin;
                raw = buf(a + 1:a + chunkBytes);
                x(:, k) = double(typecast(raw, fmt));
            end
            y = resample(x(:), up, down);
            tFirst = t(ii(1));
            j0 = round((b0 - tFirst) * fsOut);
            j1 = round((b1 - tFirst) * fsOut);
            j0 = max(j0, 0); j1 = max(j1, 0);
            j1 = min(j1, numel(y));
            pieces{s, ci}{end + 1} = y(j0 + 1:j1); %#ok<AGROW>
        end
    end
    if verbose
        fprintf('  read %7.0f to %7.0f s  (%.0f MB)\n', b0, b1, (oMax - oMin) / 1e6);
    end
end

data = struct();
for s = 1:nStore
    arrs = cell(1, numel(channels));
    for ci = 1:numel(channels)
        if isempty(pieces{s, ci})
            arrs{ci} = zeros(0, 1);
        else
            arrs{ci} = vertcat(pieces{s, ci}{:});
        end
    end
    n = min(cellfun(@numel, arrs));
    M = zeros(numel(channels), n);
    for ci = 1:numel(channels)
        M(ci, :) = arrs{ci}(1:n);
    end
    data.(stores{s}) = M;
end
out = struct('fs', fsOut, 't0', tStart, 'fsRaw', fsRaw, 'data', data);
end


function [tsqPath, tevPath] = findFiles(blockDir)
d = dir(fullfile(blockDir, '*.tsq'));
e = dir(fullfile(blockDir, '*.tev'));
assert(~isempty(d) && ~isempty(e), 'No .tsq/.tev pair in %s', blockDir);
tsqPath = fullfile(blockDir, d(1).name);
tevPath = fullfile(blockDir, e(1).name);
end


function rec = readIndex(tsqPath)
% Each .tsq record is 40 bytes: size, type, code, channel, sortcode,
% timestamp, offset, format, frequency.
fid = fopen(tsqPath, 'r');
cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>
raw = fread(fid, inf, '*uint8');
n = floor(numel(raw) / 40);
raw = reshape(raw(1:n * 40), 40, n);
rec.size      = double(typecast(reshape(raw(1:4, :), [], 1), 'int32'));
rec.type      = double(typecast(reshape(raw(5:8, :), [], 1), 'int32'));
rec.code      = double(typecast(reshape(raw(9:12, :), [], 1), 'uint32'));
rec.channel   = double(typecast(reshape(raw(13:14, :), [], 1), 'uint16'));
rec.sortcode  = double(typecast(reshape(raw(15:16, :), [], 1), 'uint16'));
rec.timestamp = typecast(reshape(raw(17:24, :), [], 1), 'double');
rec.offset    = double(typecast(reshape(raw(25:32, :), [], 1), 'int64'));
rec.format    = double(typecast(reshape(raw(33:36, :), [], 1), 'int32'));
rec.frequency = double(typecast(reshape(raw(37:40, :), [], 1), 'single'));
end


function name = formatName(code)
names = {'single', 'int32', 'int16', 'int8', 'double'};
name = names{code + 1};
end


function b = formatBytes(code)
sizes = [4 4 2 1 8];
b = sizes(code + 1);
end
