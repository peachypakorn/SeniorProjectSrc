function [synth, depthLo, depthRo] = laplLagr(img1, img2, nLev, nA, dA, dsrLagr, maxDLagr, dReg, depthLo, depthRo)
    nR = size(img1, 1);
    nC = size(img1, 2);
    nCol = size(img1, 3);
    synth = zeros([ size(img1)  nA*2]);

    if isscalar(depthLo) || isscalar(depthRo),
        [depthLo, depthRo] = lagrDepthEstimate(imresize(img1, dsrLagr), imresize(img2, dsrLagr), maxDLagr, dReg);
    end

    depthL = imresize(depthLo, [nR, nC]) / dsrLagr;
    depthR = imresize(depthRo, [nR, nC]) / dsrLagr;

    pyrL = buildLaplacianPyr(img1, nLev, 1/3);
    pyrR = buildLaplacianPyr(img2, nLev, 1/3);

    nC1 = size(pyrL{1}, 2);
    nAL = floor((nC1 - nC) / 2);
    tmp = zeros([nR, nC1]);
    tmp(:, nAL+1:nAL+nC) = depthL;
    depthL = tmp;
    tmp = zeros([nR, nC1]);
    tmp(:, nAL+1:nAL+nC) = depthR;
    depthR = tmp;

    for iA = 1:nA,
        amp = -0.5 + (iA - 0.5)*dA;

        pyrNL = {};
        for iLev = 1:nLev+1,
            cLen = size(pyrL{iLev}, 2);
            depthResize = imresize(depthL, [nR, cLen]) * cLen / nC1;
            pyrNL{iLev} = warpImage(pyrL{iLev}, - depthResize * amp, depthResize);
        end
        pyrNR = {};
        for iLev = 1:nLev+1,
            cLen = size(pyrR{iLev}, 2);
            depthResize = imresize(depthR, [nR, cLen]) * cLen / nC1;
            pyrNR{iLev} = warpImage(pyrR{iLev}, depthResize * amp, depthResize);
        end

        synth(:, :, :, nA + 1 - iA) = reconLaplacianPyr(pyrNL, nC, 1/3);
        synth(:, :, :, nA + iA) = reconLaplacianPyr(pyrNR, nC, 1/3);
    end
end

function pyr = buildLaplacianPyr(img, nLev, a)
    nR = size(img, 1);
    nC = size(img, 2);
    nCol = size(img, 3);
    nC1 = ceil((nC-1) / (2^nLev)) * (2^nLev) + 1;
    nAL = floor((nC1 - nC) / 2);
    nAR = nC1 - nC - nAL;
    tmp = zeros([nR, nC1, nCol]);
    tmp(:, nAL + 1 : nAL + nC, :) = img;
    tmp(:, 1:nAL, :) = repmat(img(:, 1, :), [1, nAL, 1]);
    tmp(:, nAL + nC + 1:nC1, :) = repmat(img(:, nC, :), [1, nAR, 1]);
    pyr{1} = tmp;
    for iLev = 1:nLev,
        cLen = ceil(nC1 / (2^iLev));
        pyr{iLev+1} = zeros([nR, cLen, nCol]);
        pyr{iLev+1}(:, 2:cLen, :) = pyr{iLev+1}(:, 2:cLen, :) + (0.25 - a/2) * pyr{iLev}(:, 1:2:cLen*2-3, :);
        pyr{iLev+1}(:, 2:cLen, :) = pyr{iLev+1}(:, 2:cLen, :) + 0.25 * pyr{iLev}(:, 2:2:cLen*2-2, :);
        pyr{iLev+1} = pyr{iLev+1} + a * pyr{iLev}(:, 1:2:cLen*2-1, :);
        pyr{iLev+1}(:, 1:cLen-1, :) = pyr{iLev+1}(:, 1:cLen-1, :) + 0.25 * pyr{iLev}(:, 2:2:cLen*2-2, :);
        pyr{iLev+1}(:, 1:cLen-1, :) = pyr{iLev+1}(:, 1:cLen-1, :) + (0.25 - a/2) * pyr{iLev}(:, 3:2:cLen*2-1, :);
        tmp = zeros(size(pyr{iLev}));
        tmp(:, 1:2:cLen*2-3, :) = tmp(:, 1:2:cLen*2-3, :) + (0.25 - a/2) * pyr{iLev+1}(:, 2:cLen, :);
        tmp(:, 2:2:cLen*2-2, :) = tmp(:, 2:2:cLen*2-2, :) + 0.25 * pyr{iLev+1}(:, 2:cLen, :);
        tmp(:, 1:2:cLen*2-1, :) = tmp(:, 1:2:cLen*2-1, :) + a * pyr{iLev+1};
        tmp(:, 2:2:cLen*2-2, :) = tmp(:, 2:2:cLen*2-2, :) + 0.25 * pyr{iLev+1}(:, 1:cLen-1, :);
        tmp(:, 3:2:cLen*2-1, :) = tmp(:, 3:2:cLen*2-1, :) + (0.25 - a/2) * pyr{iLev+1}(:, 1:cLen-1, :);
        pyr{iLev} = pyr{iLev} - tmp * 2;
    end
end

function img = reconLaplacianPyr(pyr, nC, a)
    nLev = length(pyr);
    for iLev = nLev-1:-1:1,
        cLen = size(pyr{iLev+1}, 2);
        pyr{iLev}(:, 1:2:cLen*2-3, :) = pyr{iLev}(:, 1:2:cLen*2-3, :) + 2*(0.25 - a/2) * pyr{iLev+1}(:, 2:cLen, :);
        pyr{iLev}(:, 2:2:cLen*2-2, :) = pyr{iLev}(:, 2:2:cLen*2-2, :) + 0.5*pyr{iLev+1}(:, 2:cLen, :);
        pyr{iLev}(:, 1:2:cLen*2-1, :) = pyr{iLev}(:, 1:2:cLen*2-1, :) + 2*a * pyr{iLev+1};
        pyr{iLev}(:, 2:2:cLen*2-2, :) = pyr{iLev}(:, 2:2:cLen*2-2, :) + 0.5*pyr{iLev+1}(:, 1:cLen-1, :);
        pyr{iLev}(:, 3:2:cLen*2-1, :) = pyr{iLev}(:, 3:2:cLen*2-1, :) + 2*(0.25 - a/2) * pyr{iLev+1}(:, 1:cLen-1, :);
    end
    nC1 = size(pyr{1}, 2);
    nAL = floor((nC1 - nC) / 2);
    img = pyr{1}(:, nAL+1:nAL+nC, :);
end