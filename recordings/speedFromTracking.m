function [tg, v] = speedFromTracking(tr, pxPerCm, smoothS)
%SPEEDFROMTRACKING Running speed in cm/s from the LED midpoint.
%   [TG, V] = SPEEDFROMTRACKING(TR) interpolates onto a regular 30 Hz grid and
%   smooths lightly. Frames where one tracker lost its LED read as (0,0) and
%   fall back to the other LED.

if nargin < 2 || isempty(pxPerCm), pxPerCm = 3.68; end
if nargin < 3 || isempty(smoothS), smoothS = 0.25; end
t = tr.t(:);
x = 0.5 * (tr.x1(:) + tr.x2(:));
y = 0.5 * (tr.y1(:) + tr.y2(:));
bad1 = tr.x1(:) == 0 & tr.y1(:) == 0;
bad2 = tr.x2(:) == 0 & tr.y2(:) == 0;
x(bad1) = tr.x2(bad1); y(bad1) = tr.y2(bad1);
x(bad2) = tr.x1(bad2); y(bad2) = tr.y1(bad2);
good = ~(bad1 & bad2);
tg = (t(1):1/30:t(end)).';
xi = interp1(t(good), x(good), tg, 'linear', 'extrap');
yi = interp1(t(good), y(good), tg, 'linear', 'extrap');
k = max(round(smoothS * 30), 1);
w = ones(k, 1) / k;
xi = conv(xi, w, 'same');
yi = conv(yi, w, 'same');
v = hypot(gradient(xi, tg), gradient(yi, tg)) / pxPerCm;
end
