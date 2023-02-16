function y = AddWindow(x,interval)
%AddWindow A simple function to add a Kaiser window specified by interval to x
%   x is row vector or matrix, interval is a vector with two elements,
%   marking the colomns of the beginning and ending of the window.
    w = kaiser(interval(2) - interval(1) + 1, 2.5);
    w = [zeros(1, interval(1) - 1), w'];
    w = [w, zeros(1, size(x, 2) - numel(w))];
    y = x .* w;
end

