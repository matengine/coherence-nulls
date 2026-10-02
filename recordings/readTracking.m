function tr = readTracking(blockDir)
%READTRACKING Parse the two-LED tracking log written beside a TDT tank.
%   TR = READTRACKING(BLOCKDIR) returns a struct with t, x1, y1, x2, y2.

path = fullfile(blockDir, 'tracking.txt');
txt = fileread(path);
pat = ['\(at ([0-9.]+)s\).*?T1: region:\S+ heading:\S+ x:([0-9.\-]+), y:([0-9.\-]+)\)' ...
       '.*?T2: region:\S+ heading:\S+ x:([0-9.\-]+), y:([0-9.\-]+)\)'];
tok = regexp(txt, pat, 'tokens');
n = numel(tok);
v = zeros(n, 5);
for i = 1:n
    v(i, :) = str2double(tok{i});
end
tr = struct('t', v(:, 1), 'x1', v(:, 2), 'y1', v(:, 3), 'x2', v(:, 4), 'y2', v(:, 5));
end
