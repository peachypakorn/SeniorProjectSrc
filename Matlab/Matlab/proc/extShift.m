function res = extShift(array, shift)
    nC = size(array, 2);
    res = zeros(size(array));
    if shift >= 0,
        res(:, 1+shift:nC, :) = array(:, 1:nC-shift, :);
    else,
        res(:, 1:nC+shift, :) = array(:, 1-shift:nC, :);
    end
end