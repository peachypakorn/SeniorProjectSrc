function [filteredL, filteredR] = trilateralDepthFilter(depthL, depthR, imageL, imageR, kernelR, kernelC, sigC, nItr, colProc)
    tLR = 2;
    tCol = 5 / 255.0;

    nR = size(imageL, 1);
    nC = size(imageL, 2);
    nCol = size(imageL, 3);

    depthL = round(depthL);
    depthR = round(depthR);

    hwR = floor(kernelR/2);
    hwC = floor(kernelC/2);
    %kernel = sigL*2;
    %hw = floor(kernel/2);

    %sigmapL = zeros(nR, nC);
    %sigmapR = zeros(nR, nC);
    %for iCol = 1:nCol,
    %    sigmapL(1:nR-1, 1:nC) = or(sigmapL(1:nR-1, 1:nC), ...
    %        abs(imageL(1:nR-1, 1:nC, iCol) - imageL(2:nR, 1:nC, iCol)) > tCol);
    %    sigmapL(1:nR-1, 1:nC) = or(sigmapL(1:nR-1, 1:nC), ...
    %        abs(imageL(1:nR-1, 1:nC, iCol) - imageL(2:nR, 1:nC, iCol)) > tCol);
    %    sigmapL(1:nR, 1:nC-1) = or(sigmapL(1:nR, 1:nC-1), ...
    %        abs(imageL(1:nR, 1:nC-1, iCol) - imageL(1:nR, 2:nC, iCol)) > tCol);
    %    sigmapL(1:nR, 1:nC-1) = or(sigmapL(1:nR, 1:nC-1), ...
    %        abs(imageL(1:nR, 1:nC-1, iCol) - imageL(1:nR, 2:nC, iCol)) > tCol);
    %    sigmapR(1:nR-1, 1:nC) = or(sigmapR(1:nR-1, 1:nC), ...
    %        abs(imageR(1:nR-1, 1:nC, iCol) - imageR(2:nR, 1:nC, iCol)) > tCol);
    %    sigmapR(1:nR-1, 1:nC) = or(sigmapR(1:nR-1, 1:nC), ...
    %        abs(imageR(1:nR-1, 1:nC, iCol) - imageR(2:nR, 1:nC, iCol)) > tCol);
    %    sigmapR(1:nR, 1:nC-1) = or(sigmapR(1:nR, 1:nC-1), ...
    %        abs(imageR(1:nR, 1:nC-1, iCol) - imageR(1:nR, 2:nC, iCol)) > tCol);
    %    sigmapR(1:nR, 1:nC-1) = or(sigmapR(1:nR, 1:nC-1), ...
    %        abs(imageR(1:nR, 1:nC-1, iCol) - imageR(1:nR, 2:nC, iCol)) > tCol);
    %end

    %sigmapLEx = zeros(nR+hw*2, nC+hw*2);
    %sigmapREx = zeros(nR+hw*2, nC+hw*2);
    %sigmapLEx(hw+(1:nR), hw+(1:nC)) = sigmapL;
    %sigmapREx(hw+(1:nR), hw+(1:nC)) = sigmapR;

    %sigmapLEx = sigmapLEx * (sigS - sigL) + sigL;
    %sigmapREx = sigmapREx * (sigS - sigL) + sigL;

    imageLEx = zeros(nR+hwR*2, nC+hwC*2, nCol);
    imageREx = zeros(nR+hwR*2, nC+hwC*2, nCol);
    imageLEx(hwR+(1:nR), hwC+(1:nC), :) = imageL;
    imageREx(hwR+(1:nR), hwC+(1:nC), :) = imageR;

    for itr = 1:nItr,
        [r, c] = ind2sub([nR, nC], find(ones(nR, nC)));
        c1 = c - reshape(depthL, [nR*nC, 1]);
        c1(c1<1) = 1;
        c1(c1>nC) = nC;
        checkL = abs(depthL - reshape(depthR(sub2ind([nR, nC], r, c1)), [nR, nC]));
        wL = 1 ./ checkL;
        wL(checkL > tLR) = 0;
        wL(checkL == 0) = 1;

        c2 = c + reshape(depthR, [nR*nC, 1]);
        c2(c2>nC) = nC;
        c2(c2<1) = 1;
        checkR = abs(depthR - reshape(depthL(sub2ind([nR, nC], r, c2)), [nR, nC]));
        wR = 0.5 ./ checkR;
        wR(checkR > tLR) = 0;
        wR(checkR == 0) = 1;

        %getShiftWeightL = @(dR, dC, rSt, rEd) getShiftWeight(imageLEx, sigmapLEx, sigC, nR, nC, hw, dR, dC, rSt, rEd);
        getShiftWeightL = @(dR, dC, cSt, cEd) getShiftWeightFlat(imageLEx, sigC, nR, nC, hwR, hwC, dR, dC, cSt, cEd);
        depthL = medianFilter(depthL, wL, getShiftWeightL, kernelR, kernelC, colProc);
        %imshow((8*depthL+128)/255);
        %figure;
        %getShiftWeightR = @(dR, dC, rSt, rEd) getShiftWeight(imageREx, sigmapREx, sigC, nR, nC, hw, dR, dC, rSt, rEd);
        getShiftWeightR = @(dR, dC, cSt, cEd) getShiftWeightFlat(imageREx, sigC, nR, nC, hwR, hwC, dR, dC, cSt, cEd);
        depthR = medianFilter(depthR, wR, getShiftWeightR, kernelR, kernelC, colProc);
        %imshow((8*depthR+128)/255);
        %figure;
    end

    filteredL = depthL;
    filteredR = depthR;
end

function sw = getShiftWeightFlat(imageEx, sigC, nR, nC, hwR, hwC, dR, dC, cSt, cEd)
    wC = exp(-0.5 * (sigC^-2) * sum((imageEx(hwR+(1:nR), hwC+(cSt:cEd), :) - ...
        imageEx(hwR+dR+(1:nR), hwC+dC+(cSt:cEd), :)) .^ 2, 3));
    sw = wC;
end

function sw = getShiftWeight(imageEx, sigmapEx, sigC, nR, nC, hw, dR, dC, rSt, rEd)
    wC = exp(-0.5 * (sigC^-2) * sum((imageEx(hw+(rSt:rEd), hw+(1:nC), :) - ...
        imageEx(hw+dR+(rSt:rEd), hw+dC+(1:nC), :)) .^ 2, 3));
    wD = exp(-0.5 * (sigmapEx(hw+(rSt:rEd), hw+(1:nC)) .^ -2) * (dR^2+dC^2));
    sw = wC .* wD;
end