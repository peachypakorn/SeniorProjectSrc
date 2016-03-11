function [depth] = depthUpsample(depthL, depthLPrv, image, imageL, imageLPrv, kernelR, kernelC, sigC, colProc)
    nR = size(image, 1);
    nC = size(image, 2);
    nCol = size(image, 3);
    nRL = size(imageL, 1);
    nCL = size(imageL, 2);
    hwR = floor(kernelR/2);
    hwC = floor(kernelC/2);
    nKE = (hwR*2+1) * (hwC*2+1);

    weightEx = zeros([nRL+2*hwR+2, nCL+2*hwC+2]);
    weightEx(hwR+1+(1:nRL), hwC+1+(1:nCL)) = 1;

    if length(depthLPrv) > 1,
        depthList = zeros([nRL+2*hwR+2, nCL+2*hwC+2, 2]);
        depthList(hwR+1+(1:nRL), hwC+1+(1:nCL), 1) = depthL;
        depthList(hwR+1+(1:nRL), hwC+1+(1:nCL), 2) = depthLPrv;
        imageList = zeros([nRL+2*hwR+2, nCL+2*hwC+2, nCol, 2]);
        imageList(hwR+1+(1:nRL), hwC+1+(1:nCL), :, 1) = imageL;
        imageList(hwR+1+(1:nRL), hwC+1+(1:nCL), :, 2) = imageLPrv;
    else,
        depthList = zeros([nRL+2*hwR+2, nCL+2*hwC+2, 1]);
        depthList(hwR+1+(1:nRL), hwC+1+(1:nCL), 1) = depthL;
        imageList = zeros([nRL+2*hwR+2, nCL+2*hwC+2, nCol, 1]);
        imageList(hwR+1+(1:nRL), hwC+1+(1:nCL), :, 1) = imageL;
    end

    depth = zeros(nR, nC);

    for cSt = 1:colProc:nC,
        %cSt
        cEd = min(cSt+colProc-1, nC);
        cSize = cEd-cSt+1;
        idx = (1:cSize*nR)';

        [centerR, centerC] = ind2sub([nR, cSize], idx);
        centerC = centerC + cSt - 1;
        centerR = round((centerR-0.5) * nRL / nR);
        centerC = round((centerC-0.5) * nCL / nC);

        sumD = zeros(cSize*nR, 1);
        sumW = zeros(cSize*nR, 1);
        for imgId = 1:size(depthList, 3);
            val = reshape(image(1:nR, cSt:cEd, :), [cSize*nR, nCol]);
            depthEx = depthList(:, :, imgId);
            imageLEx = imageList(:, :, :, imgId);
            for dR = -hwR:hwR,
                for dC = -hwC:hwC,
                    val2 = zeros(cSize*nR, nCol);
                    for iCol = 1:nCol,
                        val2(:, iCol) = imageLEx(sub2ind([nRL+2*hwR+2, nCL+2*hwC+2, nCol], ...
                            centerR+hwR+1+dR, centerC+hwC+1+dC, iCol * ones(cSize*nR, 1)));
                    end
                    wc = exp(-0.5 * (sigC^-2) * sum((val - val2).^2, 2));
                    ww = weightEx(sub2ind([nRL+2*hwR+2, nCL+2*hwC+2], ...
                            centerR+hwR+1+dR, centerC+hwC+1+dC));
                    dv = depthEx(sub2ind([nRL+2*hwR+2, nCL+2*hwC+2], ...
                            centerR+hwR+1+dR, centerC+hwC+1+dC));
                    weight = wc .* ww;

                    sumD = sumD + dv .* weight;
                    sumW = sumW + weight;
                end
            end
            sumD = sumD / 4;
            sumW = sumW / 4;
        end
        depth(:, cSt:cEd) = reshape(sumD ./ sumW, [nR, cSize]);
    end

    depth = depth * nC / nCL;
end