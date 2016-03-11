function res = posFreq1D(im)
    res = zeros(size(im));
    n2 = size(im, 2);
    freq = fft(im, [], 2);
    freq(:, ceil(n2/2)+1:n2, :) = 0;
    freq(:, 1, :) = 0;
    res = ifft(freq, [], 2);
end