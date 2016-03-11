function [filtered] = adaptiveGaussian(image, sigma, kernel)
    nR = size(image, 1);
    nC = size(image, 2);
    nCol = size(image, 3);

    hw = floor(kernel/2);

    imageEx = zeros(nR+hw*2, nC+hw*2, nCol);
    imageEx(hw+(1:nR), hw+(1:nC)) = image;
    sigmaEx = 0.001 * ones(nR+hw*2, nC+hw*2);
    sigmaEx(hw+(1:nR), hw+(1:nC)) = sigma;

    sumW = zeros(nR, nC, nCol);
    sumI = zeros(nR, nC, nCol);

    skip = ceil(hw / 10);
    nval = floor(hw / skip);
    for dR = -skip*nval:skip:skip*nval,
        %dR
        for dC = -skip*nval:skip:skip*nval,
            val = imageEx(hw+dR+(1:nR), hw+dC+(1:nC), :);
            sig = sigmaEx(hw+(1:nR), hw+(1:nC));
            w = repmat(exp(-0.5 * (sig .^ -2) * (dR^2+dC^2)), [1, 1, nCol]);
            sumW = sumW + w;
            sumI = sumI + w .* val;
        end
    end

    filtered = sumI ./ sumW;
end