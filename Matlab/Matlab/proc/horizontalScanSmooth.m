function [filtered] = horizontalScanSmooth(depth, moveFactor, maxHole)
    nR = size(depth, 1);
    nC = size(depth, 2);
    filtered = zeros(nR, nC);

    depthGauss = imfilter(depth, fspecial('gaussian', [5, 5], 1), 'symmetric');

    if moveFactor < 0
        cRange = nC:-1:1;
    else,
        cRange = 1:nC;
    end
    maxStep = maxHole / max(1e-6, abs(moveFactor));

    nC = 0;
    alpha = 1.0/6;
    curMin = 1e6*ones(nR, 1);
    for iC = cRange,
        if nC * alpha < 1,
            updMin = (curMin * nC + depthGauss(:, iC)) / (nC+1);
            nC = nC + 1;
        else,
            updMin = curMin * (1-alpha) + depthGauss(:, iC) * alpha;
        end
        filtered(:, iC) = min(depthGauss(:, iC), curMin + maxStep);
        curMin = min(curMin + maxStep, updMin);
    end
end