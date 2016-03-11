function [synth, depthLo, depthRo] = lagrangian(img1, img2, depthLo, depthRo, testParam, runParam)
    rho = runParam.rho;
    nA = runParam.nA;
    dA = runParam.dA;
    dsrLagr = runParam.lagrResize;
    maxDLagr = ceil(testParam.maxDepth * dsrLagr);
    if isfield(runParam, 'dReg'),
        dReg = runParam.dReg;
    else,
        dReg = 1e10;
    end
    if isfield(runParam, 'remap'),
        maxDepth = runParam.remapDepth;
    else
        maxDepth = 0;
    end
    nPP = runParam.nPushPull;

    nR = size(img1, 1);
    nC = size(img1, 2);
    nCol = size(img1, 3);
    synth = zeros([ size(img1)  nA*2]);

    if isscalar(depthLo) || isscalar(depthRo),
        [depthLo, depthRo] = lagrDepthEstimate(imresize(img1, dsrLagr), imresize(img2, dsrLagr), maxDLagr, dReg);
    end
    depthL = imresize(depthLo, [nR, nC]) / dsrLagr;
    depthR = imresize(depthRo, [nR, nC]) / dsrLagr;

    if maxDepth > 0,
        depthClp1 = maxDepth * atan(depthL / maxDepth);
        depthClp2 = maxDepth * atan(depthR / maxDepth);
        %depthClp = clip(depth, -maxDepth, maxDepth);
    else
        depthClp1 = depthL;
        depthClp2 = depthR;
    end

    if false,
        for itr = 1:30,
            h = fspecial('gaussian', [7, 7], 3);
            depthL = bsxfun(@min, depthL, imfilter(depthL, h));
            depthR = bsxfun(@min, depthR, imfilter(depthR, h));
        end
    end

    if rho == 0,
        for iA = nA:-1:1,
            movL = -depthClp1 * (iA - 0.5)*dA + depthL * 0.5;
            movR = depthClp2 * (iA - 0.5)*dA - depthR * 0.5;
            synth(:, :, :, nA + 1 - iA) = warpImage(img1, movL, depthL, nPP);
            synth(:, :, :, nA + iA) = warpImage(img2, movR, depthR, nPP);
        end
    else
        nAd = ceil((nA - 0.5) * 3 + 2 * rho * 3 + 0.5);
        synthd = zeros([ size(img1)  nAd*2]);
        for iA = nAd:-1:1,
            movL = -depthClp1 * (iA - 0.5)*dA/3 + depthL * 0.5;
            movR = depthClp2 * (iA - 0.5)*dA/3 - depthR * 0.5;
            synthd(:, :, :, nAd + 1 - iA) = warpImage(img1, movL, depthL, nPP);
            synthd(:, :, :, nAd + iA) = warpImage(img2, movR, depthR, nPP);
        end
        for iA = 1:2*nA,
            center = (iA - 1)* 3 + 1 + (2*nAd - (2*nA-1) * 3-1) / 2;
            hw = ceil(2 * rho * 3);
            w0 = 0;
            for qA = center-hw:center+hw
                dw = exp(-0.5*(qA - center)^2/(rho * 3)^2);
                synth(:, :, :, iA) = synth(:, :, :, iA) + synthd(:, :, :, qA) * dw;
                w0 = w0 + dw;
            end
            synth(:, :, :, iA) = synth(:, :, :, iA) / w0;
        end
    end


end