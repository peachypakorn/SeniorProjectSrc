function [synth, mem, phsCmpO] = eulaOnLagrComp(img1, img2, depthLo, depthRo, testParam, runParam, mem, iA, noEula, phsCmpO)
    rho = runParam.rho;
    ht = min(runParam.ht, ceil(log(size(img1, 2)) / log(2)) - 4);
    nA = runParam.nA;
    dA = runParam.dA;
    dsrLagr = runParam.eurLagrResize;
    vis = struct();
    if isfield(runParam, 'dReg'),
        dReg = runParam.dReg;
    else,
        dReg = 1e10;
    end
    if isfield(runParam, 'forcematch') && runParam.forcematch,
        sigWrap = runParam.sigWrap;
        appWrap = runParam.appWrap;
    else,
        sigWrap = 0;
        appWrap = 0;
    end
    if isfield(runParam, 'recR'),
        recR = runParam.recR;
    else,
        recR = 0;
    end
    if isfield(runParam, 'expSinc'),
        expSinc = runParam.expSinc;
    else,
        expSinc = false;
    end
    if isfield(runParam, 'expShear'),
        expShear = runParam.expShear;
    else,
        phsCmpO = 0;
        expShear = false;
    end
    if isfield(runParam, 'expBlend'),
        expBlend = runParam.expBlend;
        blendLayer = runParam.blendLayer;
    else,
        blendLayer = 0;
        expBlend = false;
    end
    if isfield(runParam, 'zipzapBlend'),
        zipzapBlend = runParam.zipzapBlend;
        blendPeriod = runParam.blendPeriod;
        blendAmp = runParam.blendAmp;
    else,
        blendPeriod = 0;
        zipzapBlend = false;
        blendAmp = 0;
    end
    if isfield(runParam, 'expPrealign'),
        expPrealign = runParam.expPrealign;
    else,
        expPrealign = false;
    end
    if isfield(runParam, 'avgPrealign'),
        avgPrealign = runParam.avgPrealign;
    else,
        avgPrealign = false;
    end
    if isfield(runParam, 'twoSideInterp'),
        twoSideInterp = runParam.twoSideInterp;
    else,
        twoSideInterp = false;
    end
    if isfield(runParam, 'remap'),
        remap = runParam.remap;
        if remap,
            origMinD = runParam.origMinD;
            origMaxD = runParam.origMaxD;
            gama = runParam.gama;
            fixZero = runParam.fixZero;
        end
    else,
        remap = false;
    end

    waveletExt = 4;

    nR = size(img1, 1);
    nC = size(img1, 2);
    nCol = size(img1, 3);
    maxAmp = (nA-1)*dA;
    sigma_v_fac = 1e-6;

    lagrDepthL = depthLo;
    lagrDepthR = depthRo;
    dsrLagr = size(lagrDepthL, 2) / size(img1, 2);

    %h = fspecial('gaussian', [5, 5], 1);
    %lagrDepthL = imfilter(lagrDepthL, h, 'symmetric');
    %lagrDepthR = imfilter(lagrDepthR, h, 'symmetric');


    sizFilt = 2.^[0:-1:-ht];
    filters = get1DRadialFilters(nC, sizFilt, 1);
    nLev = max(size(filters));
    depPhsRatios = getDepthPhaseRatio(filters);

    pind = buildBandIndices(filters, 1);
    nF = pind(nLev, 1) + pind(nLev, 2) - 1;

    pindRef = buildBandIndices(filters, waveletExt);
    nFRef = pindRef(nLev, 1) + pindRef(nLev, 2) - 1;
    if ~isstruct(mem),
        pyr1 = zeros([nF nR nCol]);
        pyr2 = zeros([nF nR nCol]);

        refInd1 = ones([nF nR nCol]);
        refInd2 = ones([nF nR nCol]);
        preDepth1 = zeros([nF nR]);
        preDepth2 = zeros([nF nR]);
        preDepthCont1 = zeros([nF nR]);
        preDepthCont2 = zeros([nF nR]);
        preDepthRef1 = zeros([nFRef nR]);
        preDepthRef2 = zeros([nFRef nR]);

        for iR=1:nR,
            for iCol = 1:nCol,
                pyr1(:, iR, iCol) = buildSCFpyrGen1D(img1(iR, :, iCol), filters, 1);
                pyr2(:, iR, iCol) = buildSCFpyrGen1D(img2(iR, :, iCol), filters, 1);
            end
        end;

        for iLev = 1:nLev-1,
            ind = pyrBandIndices(pind, iLev);
            indRef = pyrBandIndices(pindRef, iLev);
            nO = pind(iLev, 2);
            nRef = pindRef(iLev, 2);
            if nO == nRef,
                refStart = 0.5;
            else,
                refStart = 0;
            end

            if ~expPrealign && ~avgPrealign,
                levDepthL = imresize(lagrDepthL, [nR, pind(iLev, 2)*4])'/dsrLagr;
                sigma_v = sigma_v_fac;
                levDepthL = reduce(@min, levDepthL, 4, 0);
                preDepthCont1(ind, :) = levDepthL;

                levDepthR = imresize(lagrDepthR, [nR, pind(iLev, 2)*4])'/dsrLagr;
                levDepthR = reduce(@min, levDepthR, 4, 0);
                preDepthCont2(ind, :) = levDepthR;
            elseif expPrealign,
                p = -nC / 40;
                lagrDepthL = imresize(lagrDepthL, [nR, nC])/dsrLagr;
                levDepthL = exp(lagrDepthL'/p);
                levDepthL = p * log(reduceToAvg(levDepthL, [pind(iLev, 2), nR], 2));
                preDepthCont1(ind, :) = levDepthL;

                lagrDepthR = imresize(lagrDepthR, [nR, nC])/dsrLagr;
                levDepthR = exp(lagrDepthR'/p);
                levDepthR = p * log(reduceToAvg(levDepthR, [pind(iLev, 2), nR], 2));
                preDepthCont2(ind, :) = levDepthR;
                imshow(imresize(levDepthL', [nR, nC]) / 80 + 0.5);
                figure
            else
                levDepthL = imresize(lagrDepthL, [nR, nC])'/dsrLagr;
                levDepthL = reduceToAvg(levDepthL, [pind(iLev, 2), nR], 1);
                levDepthRefL = imresize(lagrDepthL, [nR, nC])'/dsrLagr;
                levDepthRefL = reduceToAvg(levDepthRefL, [pindRef(iLev, 2), nR], pindRef(iLev, 2) / pind(iLev, 2));
                preDepthCont1(ind, :) = levDepthL;
                preDepthRef1(indRef, :) = levDepthRefL;

                levDepthR = imresize(lagrDepthR, [nR, nC])'/dsrLagr;
                levDepthR = reduceToAvg(levDepthR, [pind(iLev, 2), nR], 1);
                levDepthRefR = imresize(lagrDepthR, [nR, nC])'/dsrLagr;
                levDepthRefR = reduceToAvg(levDepthRefR, [pindRef(iLev, 2), nR], pindRef(iLev, 2) / pind(iLev, 2));
                preDepthCont2(ind, :) = levDepthR;
                preDepthRef2(indRef, :) = levDepthRefR;
            end
            levDepthL(isnan(levDepthL)) = 0;
            preDepth1(ind, :) = repmat((-0.5 + (1:nO)')*(nC/nO), [1, nR]) - levDepthL;
            preDepth1(ind, :) = clip(round(preDepth1(ind, :) * (nRef / nC) - refStart), 0, nRef-1);
            refInd1(ind, :, :) = repmat(preDepth1(ind, :) + pindRef(iLev, 1), [1, 1, nCol]);
            preDepth1(ind, :) = (preDepth1(ind, :) + refStart) * (nC / nRef) - ...
                    repmat((-0.5 + (1:nO)')*(nC/nO), [1, nR]);

            levDepthR(isnan(levDepthR)) = 0;
            preDepth2(ind, :) = repmat((-0.5 + (1:nO)')*(nC/nO), [1, nR]) + levDepthR;
            preDepth2(ind, :) = clip(round(preDepth2(ind, :) * (nRef / nC) - refStart), 0, nRef-1);
            refInd2(ind, :, :) = repmat(preDepth2(ind, :) + pindRef(iLev, 1), [1, 1, nCol]);
            preDepth2(ind, :) = repmat((-0.5 + (1:nO)')*(nC/nO), [1, nR]) - ...
                    (preDepth2(ind, :) + refStart) * (nC / nRef);
        end
        %% Computing Filter
        phsDif1 = zeros([nF nR nCol]);
        phsDif2 = zeros([nF nR nCol]);
        weight1 = zeros([nF nR nCol]);
        weight2 = zeros([nF nR nCol]);
        pyr1R = zeros([nFRef nR nCol]);
        pyr2R = zeros([nFRef nR nCol]);

        for iR=1:nR,
            for iCol = 1:nCol,
                pyr1Ref = buildSCFpyrGen1D(img1(iR, :, iCol), filters, waveletExt);
                pyr2Ref = buildSCFpyrGen1D(img2(iR, :, iCol), filters, waveletExt);
                pyr1R(:, iR, iCol) = pyr1Ref;
                pyr2R(:, iR, iCol) = pyr2Ref;
            end
        end

        pyr1Rem = pyr1;
        pyr2Rem = pyr2;
        if 0
            for iLev = nLev:-1:1,
                ind = pyrBandIndices(pind, iLev);
                sigma_v = 2^(iLev-2);
                kV = ceil(4*sigma_v);
                h_v = fspecial('gaussian', [1, kV], sigma_v);
                pyr1(ind, :, :) = imfilter(pyr1(ind, :, :), h_v, 'symmetric');
                pyr2(ind, :, :) = imfilter(pyr2(ind, :, :), h_v, 'symmetric');
                pyr1R(ind, :, :) = imfilter(pyr1R(ind, :, :), h_v, 'symmetric');
                pyr2R(ind, :, :) = imfilter(pyr2R(ind, :, :), h_v, 'symmetric');
            end
        end

        depDif1 = zeros(nF, nR);
        depDif2 = zeros(nF, nR);
        for iR=1:nR,
            for iCol = 1:nCol,
                phsDif1(:, iR, iCol) = angle(pyr2R(refInd1(:, iR, iCol), iR, iCol)) - angle(pyr1(:, iR, iCol));
                phsDif2(:, iR, iCol) = angle(pyr2(:, iR, iCol)) - angle(pyr1R(refInd2(:, iR, iCol), iR, iCol));
                %weight1(:, iR, iCol) = bsxfun(@min, abs(pyr2Ref(refInd1(:, iR, iCol))), abs(pyr1(:, iR, iCol)));
                %weight2(:, iR, iCol) = bsxfun(@min, abs(pyr2(:, iR, iCol)), abs(pyr1Ref(refInd2(:, iR, iCol))));
                weight1(:, iR, iCol) = abs(pyr1(:, iR, iCol));
                weight2(:, iR, iCol) = abs(pyr2(:, iR, iCol));
            end
            depDif1(:, iR) = abs(preDepthCont1(:, iR) - preDepthRef2(refInd1(:, iR, 1), iR));
            depDif2(:, iR) = abs(preDepthCont2(:, iR) - preDepthRef1(refInd2(:, iR, 1), iR));
        end

        %pyr1 = pyr1Rem;
        %pyr2 = pyr2Rem;

        sigma_v = sigma_v_fac;
        ind = pyrBandIndices(pind, nLev);
        h = fspecial('gaussian', [ceil(sigma_v*2)*2+1, 1], sigma_v+1e-10);
        tmpDepthL = -imfilter(imresize(lagrDepthL, [nR, pind(nLev, 2)]), h, 'symmetric')'/dsrLagr;
        preDepth1(ind, :) =  tmpDepthL;%reduce(@max, tmpDepthL, 4, 0);
        tmpDepthR = -imfilter(imresize(lagrDepthR, [nR, pind(nLev, 2)]), h, 'symmetric')'/dsrLagr;
        preDepth2(ind, :) =  tmpDepthR;%reduce(@max, tmpDepthR, 4, 0);
        wPhs1 = ones([nF, nR]);
        wPhs2 = ones([nF, nR]);
        sumPhs1 = zeros([nF, nR]);
        sumPhs2 = zeros([nF, nR]);

        expV1 = exp(phsDif1*1i) .* weight1(:, :, :);
        expV2 = exp(phsDif2*1i) .* weight2(:, :, :);
        phsDif1 = angle(sum(expV1, 3));
        phsDif2 = angle(sum(expV2, 3));
        weight1 = abs(sum(expV1, 3));
        weight2 = abs(sum(expV2, 3));
        att_c = pi / 4;
        weight1 = weight1 .* exp(-(phsDif1.^2) / (2*att_c.^2));
        weight2 = weight2 .* exp(-(phsDif2.^2) / (2*att_c.^2));
        if 0,
            for iLev = nLev-1:-1:1,
                ind = pyrBandIndices(pind, iLev);
                att_d = (pi / 64) * depPhsRatios{iLev};
                weight1(ind, :) = weight1(ind, :) .* exp(-depDif1(ind, :).^2/(2*att_d^2)) * (size(ind, 2) / nC);
                weight2(ind, :) = weight2(ind, :) .* exp(-depDif2(ind, :).^2/(2*att_d^2)) * (size(ind, 2) / nC);
            end
        end

        ind = pyrBandIndices(pind, nLev);
        phsDif1(ind, :) = 0;
        phsDif2(ind, :) = 0;
        weight1(ind, :) = 1;
        weight2(ind, :) = 1;

        for k = nLev-1:-1:1,
            ind = pyrBandIndices(pind,k);
            phsDif1(ind, :) = -preDepth1(ind, :) / depPhsRatios{k} + phsDif1(ind, :);
            phsDif2(ind, :) = -preDepth2(ind, :) / depPhsRatios{k} + phsDif2(ind, :);
            %imshow(imresize(-preDepth1(ind, :)' / 100 + 0.5, [nR, nC]));
            %figure;
            %imshow(imresize(phsDif1(ind, :)'*depPhsRatios{k} / 100 + 0.5, [nR, nC]));
            %figure;
        end

        for k = nLev-1:-1:1,
            ind = pyrBandIndices(pind,k);
            %sigma_h = max(1, sqrt(4 / (nC / length(ind))));
            sigma_h = 2; %XXX
            sigma_v = sigma_h * nC / length(ind);
            %sigma_h = 2; %XXX
            %sigma_v = 2; %XXX
            %sigma_h = sigma_v * length(ind) / nR;
            sigma_c = (pi) * depPhsRatios{k};
            kH = ceil(sigma_h*2);
            kV = ceil(sigma_v*2);
            %kH = 5; %XXX
            %kV = 5; %XXX
            %k

            if false,
                h_h = fspecial('gaussian', [kH, 1], sigma_h);
                h_v = fspecial('gaussian', [1, kV], sigma_v);
                sumPhs1(ind, :) = sum(imfilter(imfilter(phsDif1(ind, :) .* weight1(ind, :), h_h, 'circular'), h_v, 'symmetric'), 3);
                sumPhs2(ind, :) = sum(imfilter(imfilter(phsDif2(ind, :) .* weight2(ind, :), h_h, 'circular'), h_v, 'symmetric'), 3);
                wPhs1(ind, :) = sum(imfilter(imfilter(weight1(ind, :), h_h, 'circular'), h_v, 'symmetric'), 3)+1e-20;
                wPhs2(ind, :) = sum(imfilter(imfilter(weight2(ind, :), h_h, 'circular'), h_v, 'symmetric'), 3)+1e-20;
            else
                if 0,
                    sumPhs1(ind, :) = bilateral(phsDif1(ind, :), ...
                            phsDif1(ind, :), weight1(ind, :), sigma_v, sigma_h, sigma_c);
                    sumPhs2(ind, :) = bilateral(phsDif2(ind, :), ...
                            phsDif2(ind, :), weight2(ind, :), sigma_v, sigma_h, sigma_c);
                elseif 1,
                    %shiftFcn = @(dR, dC, cSt, cEd) getShiftWeight(dR, dC, cSt, cEd, length(ind), sigma_h, sigma_v);
                    shiftFcn = @(dR, dC, cSt, cEd) getShiftWeightConst(dR, dC, cSt, cEd, length(ind));
                    sumPhs1(ind, :) = medianFilter(phsDif1(ind, :), weight1(ind, :), shiftFcn, kH, kV, 20);
                    sumPhs2(ind, :) = medianFilter(phsDif2(ind, :), weight2(ind, :), shiftFcn, kH, kV, 20);
                else,
                    sumPhs1(ind, :) = (boxfilter(phsDif1(ind, :) .* weight1(ind, :), [kH, kV]) + ...
                                       5e-2 * preDepthCont1(ind, :)/depPhsRatios{k}) ./ ...
                                      (boxfilter(weight1(ind, :), [kH, kV]) + 5e-2);
                    sumPhs2(ind, :) = (boxfilter(phsDif2(ind, :) .* weight2(ind, :), [kH, kV]) + ...
                                       5e-2 * preDepthCont2(ind, :)/depPhsRatios{k}) ./ ...
                                      (boxfilter(weight2(ind, :), [kH, kV]) + 5e-2);
                end
            end
            if 0%k >= nLev - 6,
                global RESULTPATH;
                pat = sprintf('%sdepthvis/', RESULTPATH);
                chkImg = zeros(nR, nC, 3);
                checkImg(:, :, 2) = imresize((sumPhs1(ind, :) ./ wPhs1(ind, :) - phsDif1(ind, :)) * depPhsRatios{k} / (2*testParam.maxDepth) + 0.5, [nC, nR])';
                checkImg(:, :, 3) = checkImg(:, :, 2);
                %checkImg(:, :, 1) = img1(:, :, 2);
                imwrite(checkImg, sprintf('%s%02d_filter_corr.jpg', pat, k), 'jpg');
                checkImg(:, :, 2) = imresize(abs(weight1(ind, :)')*10, [nR, nC]);
                checkImg(:, :, 3) = checkImg(:, :, 2);
                imwrite(checkImg, sprintf('%s%02d_weight.jpg', pat, k), 'jpg');
                checkImg(:, :, 2) = imresize(phsDif1(ind, :)' * depPhsRatios{k} / (2*testParam.maxDepth)+0.5, [nR, nC]);
                checkImg(:, :, 3) = checkImg(:, :, 2);
                imwrite(checkImg, sprintf('%s%02d_before.jpg', pat, k), 'jpg');
                checkImg(:, :, 2) = imresize(sumPhs1(ind, :)' ./ wPhs2(ind, :)' * depPhsRatios{k} / (2*testParam.maxDepth) +0.5, [nR, nC]);
                checkImg(:, :, 3) = checkImg(:, :, 2);
                imwrite(checkImg, sprintf('%s%02d_after.jpg', pat, k), 'jpg');
                checkImg(:, :, 2) = imresize(preDepth1(ind, :)'/  (2*testParam.maxDepth) +0.5, [nR, nC]);
                checkImg(:, :, 3) = checkImg(:, :, 2);
                imwrite(checkImg, sprintf('%s%02d_prealign.jpg', pat, k), 'jpg');
                %imshow(checkImg);
                %figure;
                %pause
            end
        end
        %wPhs1(:) = 1;
        %wPhs2(:) = 1;
        phsDif1 = sumPhs1 ./ wPhs1;
        phsDif2 = sumPhs2 ./ wPhs2;
        for k = nLev-1:-1:1,
            ind = pyrBandIndices(pind,k);
            %imshow(imresize(phsDif1(ind, :)'*depPhsRatios{k} / 100 + 0.5, [nR, nC]));
            %figure;
            phsDif1(ind, :) = preDepth1(ind, :)/depPhsRatios{k} + phsDif1(ind, :);
            phsDif2(ind, :) = preDepth2(ind, :)/depPhsRatios{k} + phsDif2(ind, :);

        end

        depth1 = zeros([nF, nR]);
        depth2 = zeros([nF, nR]);
        if ~noEula,
            for k = nLev:-1:1,
                ind = pyrBandIndices(pind,k);
                depth1(ind, :) = - preDepth1(ind, :) + phsDif1(ind, :) * depPhsRatios{k};
                depth2(ind, :) = - preDepth2(ind, :) + phsDif2(ind, :) * depPhsRatios{k};
            end
        else,
            depth1 = preDepthCont1;
            depth2 = preDepthCont2;
        end

        depth1(isnan(depth1)) = 0;
        depth2(isnan(depth2)) = 0;
        if remap,
            depthClp1 = remapping(depth1, remap, origMinD, origMaxD, gama, fixZero);
            depthClp2 = remapping(depth2, remap, origMinD, origMaxD, gama, fixZero);
            preDepthClp1 = remapping(preDepthCont1, remap, origMinD, origMaxD, gama, fixZero);
            preDepthClp2 = remapping(preDepthCont2, remap, origMinD, origMaxD, gama, fixZero);

            if 0,
                for k = nLev:-1:1,
                    ind = pyrBandIndices(pind,k);
                    imshow(imresize(depth1(ind, :)' / origMaxD, [nR, nC]));
                    figure;
                    imshow(imresize(depthClp1(ind, :)' / origMaxD, [nR, nC]));
                    figure;
                end
            end
        else
            depthClp1 = depth1;
            depthClp2 = depth2;
            preDepthClp1 = preDepthCont1;
            preDepthClp2 = preDepthCont2;
        end

        if expShear,
            shear1 = zeros(nC, nR);
            shear2 = zeros(nC, nR);
            shearW1 = zeros(nC, nR);
            shearW2 = zeros(nC, nR);
            for iLev = nLev-1:-1:1,
                ind = pyrBandIndices(pind, iLev);
                indRef = pyrBandIndices(pindRef, iLev);
                nO = pind(iLev, 2);
                nRef = pindRef(iLev, 2);
                if nO == nRef,
                    refStart = 0.5;
                else,
                    refStart = 0;
                end

                wDisp1 = abs(depthClp1(ind, :) * 2 * nA * dA) / depPhsRatios{iLev};
                if iLev > 1,
                    wDisp1 = wDisp1 - pi;
                else,
                    wDisp1 = wDisp1 + pi;
                end
                mask = wDisp1 > 4*pi;
                wDisp1(mask) = 8 * pi - wDisp1(mask);
                depPhsRatios{iLev}
                wDisp1 = max(wDisp1, 1e-3/depPhsRatios{iLev});
                preDepth1(ind, :) = repmat((-0.5 + (1:nO)')*(nC/nO), [1, nR]) + ...
                                    depthClp1(ind, :) * 2 * nA * dA - depth1(ind, :);
                refInd1 = clip(round(preDepth1(ind, :) * (nRef / nC) - refStart), 0, nRef-1) + pindRef(iLev, 1);
                preDepth1(ind, :) = (refInd1 - pindRef(iLev, 1) + refStart) * (nC / nRef) - preDepth1(ind, :);

                wDisp2 = abs(depthClp2(ind, :) * 2 * nA * dA) / depPhsRatios{iLev};
                if iLev > 1,
                    wDisp2 = wDisp2 - pi;
                else,
                    wDisp2 = wDisp2 + pi;
                end
                mask = wDisp2 > 4*pi;
                wDisp2(mask) = 8 * pi - wDisp2(mask);
                wDisp2 = max(wDisp2, 1e-3/depPhsRatios{iLev});
                preDepth2(ind, :) = repmat((-0.5 + (1:nO)')*(nC/nO), [1, nR]) - ...
                                    depthClp2(ind, :) * 2 * nA * dA + depth2(ind, :);
                refInd2 = clip(round(preDepth2(ind, :) * (nRef / nC) - refStart), 0, nRef-1) + pindRef(iLev, 1);
                preDepth2(ind, :) = preDepth2(ind, :) - (refInd2 - pindRef(iLev, 1) + refStart) * (nC / nRef);


                phsDif1 = zeros(nO, nR, nCol);
                phsDif2 = zeros(nO, nR, nCol);
                weight1 = zeros(nO, nR, nCol);
                weight2 = zeros(nO, nR, nCol);

                for iR=1:nR,
                    for iCol = 1:nCol,
                        phsDif1(:, iR, iCol) = angle(pyr2R(refInd1(:, iR), iR, iCol)) - angle(pyr1(ind, iR, iCol));
                        phsDif2(:, iR, iCol) = angle(pyr2(ind, iR, iCol)) - angle(pyr1R(refInd2(:, iR), iR, iCol));
                        weight1(:, iR, iCol) = bsxfun(@min, abs(pyr2R(refInd1(:, iR), iR, iCol)), abs(pyr1(ind, iR, iCol)));
                        weight2(:, iR, iCol) = bsxfun(@min, abs(pyr2(ind, iR, iCol)), abs(pyr1R(refInd2(:, iR), iR, iCol)));
                    end
                end

                phsDif1 = angle(sum(exp(phsDif1*1i) .* weight1(:, :, :), 3));
                phsDif2 = angle(sum(exp(phsDif2*1i) .* weight2(:, :, :), 3));
                weight1 = sum(weight1, 3) .* wDisp1;
                weight2 = sum(weight2, 3) .* wDisp2;

                phsDif1 = -preDepth1(ind, :) + phsDif1 * depPhsRatios{iLev};
                phsDif2 = -preDepth2(ind, :) + phsDif2 * depPhsRatios{iLev};

                %imshow(imresize(phsDif1'/16+0.5, [nR, nC]));
                %figure;

                phsDif1 = exp(i * phsDif1 / depPhsRatios{iLev});
                phsDif2 = exp(i * phsDif2 / depPhsRatios{iLev});

                sigma_h = length(ind) / 80;
                sigma_v = sigma_h * nR / length(ind);
                kH = ceil(sigma_h*4);
                kV = ceil(sigma_v*4);

                h_h = fspecial('gaussian', [kH, 1], sigma_h);
                h_v = fspecial('gaussian', [1, kV], sigma_v);
                sumPhs1(ind, :) = sum(imfilter(imfilter(phsDif1 .* weight1, h_h, 'circular'), h_v, 'symmetric'), 3);
                sumPhs2(ind, :) = sum(imfilter(imfilter(phsDif2 .* weight2, h_h, 'circular'), h_v, 'symmetric'), 3);
                wPhs1(ind, :) = sum(imfilter(imfilter(weight1, h_h, 'circular'), h_v, 'symmetric'), 3)+1e-20;
                wPhs2(ind, :) = sum(imfilter(imfilter(weight2, h_h, 'circular'), h_v, 'symmetric'), 3)+1e-20;

                phsDif1 = sumPhs1(ind, :) ./ wPhs1(ind, :);
                phsDif2 = sumPhs2(ind, :) ./ wPhs2(ind, :);
                phsDif1 = angle(phsDif1) * depPhsRatios{iLev};
                phsDif2 = angle(phsDif2) * depPhsRatios{iLev};

                imshow(imresize((depthClp1(ind, :)-phsDif1/(2 * nA * dA))'/40+0.5, [nR, nC]));
                figure;
                imshow(imresize(phsDif1'/16+0.5, [nR, nC]));
                figure;

                shear1 = shear1 + imresize(phsDif1 .* wDisp1, [nC, nR]);
                shear2 = shear2 + imresize(phsDif2 .* wDisp2, [nC, nR]);
                shearW1 = shearW1 + imresize(wDisp1, [nC, nR]);
                shearW2 = shearW2 + imresize(wDisp2, [nC, nR]);
                imshow(imresize(wDisp1', [nR, nC])/(4*pi));
                figure
            end

            shear1 = shear1 ./ shearW1;
            shear2 = shear2 ./ shearW2;
            imshow(imresize(shear1'/40+0.5, [nR, nC]));
            figure;
            imshow(imresize(shearW1'/(8*pi), [nR, nC]));
            figure;

            for k = nLev:-1:1,
                ind = pyrBandIndices(pind, k);
                nO = length(ind);
                depthClp1(ind, :) = depthClp1(ind, :) - imresize(shear1, [nO, nR])/(2 * nA * dA);
                depthClp2(ind, :) = depthClp2(ind, :) - imresize(shear2, [nO, nR])/(2 * nA * dA);
            end
        end

        phsDifAA1 = zeros([nF, nR]);
        phsDifAA2 = zeros([nF, nR]);
        for k = nLev-1:-1:1,
            ind = pyrBandIndices(pind, k);
            %imshow(imresize(depthClp1(ind, :)', [nR, nC]) / 88 + 0.5);
            %figure;
            %imshow(imresize(phsDif1(ind, :)', [nR, nC]) / (2*pi) + 0.5);
            %figure;
            phsDifAA1(ind, :) = dA * depthClp1(ind, :) / depPhsRatios{k};
            phsDifAA2(ind, :) = dA * depthClp2(ind, :) / depPhsRatios{k};
        end

        depthClp1(isnan(depthClp1)) = 0;
        depthClp2(isnan(depthClp2)) = 0;

        % Saving values
        mem = struct();
        mem.depthClp1 = depthClp1;
        mem.depthClp2 = depthClp2;
        mem.depth1 = depth1;
        mem.depth2 = depth2;
        mem.preDepthCont1 = preDepthCont1;
        mem.preDepthCont2 = preDepthCont2;
        mem.preDepthClp1 = preDepthClp1;
        mem.preDepthClp2 = preDepthClp2;
        mem.phsDifAA1 = phsDifAA1;
        mem.phsDifAA2 = phsDifAA2;
        mem.pyr1 = pyr1;
        mem.pyr2 = pyr2;
    end

    if recR == 0,
        recR = 1:nR;
    end

    for c = 1:nCol,
        %% Reconstructing
        pyrRec = zeros(nF, length(recR));
        pyrBld = zeros(nF, length(recR));
        if ~zipzapBlend,
            pyrRec = synthWrap(mem, c, nA, iA, dA, rho, recR, pind, filters, 1, twoSideInterp, depPhsRatios);
            if expBlend,
                if iA == 1 || iA == 2*nA,
                    blendLev = max(blendLayer - 2, 1);
                    blendDis = 0.5;
                elseif iA == 2 || iA == 2*nA-1,
                    blendLev = max(blendLayer - 1, 1);
                    blendDis = 1.5;
                else,
                    blendLev = -1;
                end
                if blendLev > 0,
                    if iA <= nA,
                        pyrBld = synthWrap(mem, c, nA, iA+2*nA, dA, rho, recR, pind, filters, blendLev, twoSideInterp, depPhsRatios);
                    else,
                        pyrBld = synthWrap(mem, c, nA, iA-2*nA, dA, rho, recR, pind, filters, blendLev, twoSideInterp, depPhsRatios);
                    end
                    for iLev = nLev:-1:blendLev,
                        sigma = 2^(iLev - blendLayer);
                        wt = 1 / (1 + exp(-2*blendDis / sigma));
                        ind = pyrBandIndices(pind, iLev);
                        pyrRec(ind, :) = pyrRec(ind, :) * wt + pyrBld(ind, :) * (1-wt);
                    end
                end
            end
        else
            wave = 4 * mod(recR / blendPeriod, 1);
            wave(wave > 1) = 2 - wave(wave > 1);
            wave(wave < -1) = -2 - wave(wave < -1);
            cutLoc = 0.5 + blendAmp * wave;
            blendLev = max(blendLayer - 2, 1);
            if iA <= nA,
                mask1 = cutLoc <= iA;
                pyrRec(:, mask1) = synthWrap(mem, c, nA, iA, dA, rho, recR(mask1), pind, filters, 1, twoSideInterp, depPhsRatios);
                pyrBld(:, mask1) = synthWrap(mem, c, nA, iA+2*nA, dA, rho, recR(mask1), pind, filters, blendLev, twoSideInterp, depPhsRatios);
                mask2 = cutLoc > iA;
                pyrRec(:, mask2) = synthWrap(mem, c, nA, iA+2*nA, dA, rho, recR(mask2), pind, filters, 1, twoSideInterp, depPhsRatios);
                pyrBld(:, mask2) = synthWrap(mem, c, nA, iA, dA, rho, recR(mask2), pind, filters, blendLev, twoSideInterp, depPhsRatios);
                blendDis = abs(iA - cutLoc);
            else,
                mask1 = cutLoc+2*nA >= iA;
                pyrRec(:, mask1) = synthWrap(mem, c, nA, iA, dA, rho, recR(mask1), pind, filters, 1, twoSideInterp, depPhsRatios);
                pyrBld(:, mask1) = synthWrap(mem, c, nA, iA-2*nA, dA, rho, recR(mask1), pind, filters, blendLev, twoSideInterp, depPhsRatios);
                mask2 = cutLoc+2*nA < iA;
                pyrRec(:, mask2) = synthWrap(mem, c, nA, iA-2*nA, dA, rho, recR(mask2), pind, filters, 1, twoSideInterp, depPhsRatios);
                pyrBld(:, mask2) = synthWrap(mem, c, nA, iA, dA, rho, recR(mask2), pind, filters, blendLev, twoSideInterp, depPhsRatios);
                blendDis = abs(iA - cutLoc - 2*nA);
            end
            for iLev = nLev:-1:blendLev,
                sigma = 2^(iLev - blendLayer);
                ind = pyrBandIndices(pind, iLev);
                wt = repmat(reshape(1 ./ (1 + exp(-2*blendDis / sigma)), [1, length(recR)]), [length(ind), 1]);
                pyrRec(ind, :) = pyrRec(ind, :) .* wt + pyrBld(ind, :) .* (1-wt);
            end
        end

        %else,
        %    [synth1, weight1] = ...
        %            nufftWarp(mem.pyr1(:, recR, c) .* AA1, ...
        %            -depChg12(:, recR), -mem.depthClp1(:, recR), pind, filters, 2, 4, depPhsRatios, dir1);
        %    [synth2, weight2] = ...
        %            nufftWarp(mem.pyr2(:, recR, c) .* AA2, ...
        %            -depChg22(:, recR), -mem.depthClp2(:, recR), pind, filters, 2, 4, depPhsRatios, dir2);
        %    pyrRec(:, :) = (synth1 .* weight1 + synth2 .* weight2) ./ (weight1 + weight2);
        %end

        if false%c == 1 && iA == 8,
            for iLev = nLev-1:-1:nLev-5,
                for iI = 1:length(recR),
                    cutPyr = zeros(size(pyrRec(:, iI)));
                    cutLoc = pind(iLev, 1);
                    cutPyr(cutLoc:size(cutPyr,1), :) = pyrRec(cutLoc:size(cutPyr,1), iI);
                    recons = reconSCFpyrGen1D(cutPyr, pind, filters, 1, nLev) ;
                    synth(iI, :, c) = recons;
                end
                imshow(synth(:, :, c));
                figure;
            end
        end
        pyrRec = pyrRec;

        for iI = 1:length(recR),
            recons = reconSCFpyrGen1D(pyrRec(:, iI), pind, filters, 1, nLev) ;
            synth(iI, :, c) = recons;
        end
        %pause
        %close all;
    end
end

function pyrRec = synthWrap(mem, c, nA, iA, dA, rho, recR, pind, filters, minLev, twoSideInterp, depPhsRatios)
    amp = (nA+0.5-iA)*dA;

    depChg1 = mem.depthClp1 * amp - mem.depth1/2;
    depChg2 = mem.depthClp2 * amp + mem.depth2/2;
    preDepChg1 = mem.preDepthClp1 * amp - mem.preDepthCont1/2;
    preDepChg2 = mem.preDepthClp2 * amp + mem.preDepthCont2/2;

    if amp > 0.5,
        dir1 = -1;
        dir2 = -1;
    elseif amp > -0.5,
        dir1 = 1;
        dir2 = -1;
    else,
        dir1 = 1;
        dir2 = 1;
    end

    AA1 = exp( -( mem.phsDifAA1(:, recR) * rho).^2 / 2);
    AA2 = exp( -( mem.phsDifAA2(:, recR) * rho).^2 / 2);

    %if iA <= nA,
    if (amp > 0.5 && twoSideInterp) || (amp > 0 && ~twoSideInterp),
        pyr1 = occulsion(mem.pyr1(:, recR, c), depChg1(:, recR), mem.depth1(:, recR), preDepChg1(:, recR), mem.preDepthCont1(:, recR), pind);
        pyrRec = nufftWarp(pyr1(:, :) .* AA1, ...
                depChg1(:, recR), -mem.depthClp1(:, recR), pind, filters, 2, 4, depPhsRatios, -1, minLev);
    elseif amp > -0.5 && twoSideInterp,
        pyr1 = occulsion(mem.pyr1(:, recR, c), depChg1(:, recR), mem.depth1(:, recR), preDepChg1(:, recR), mem.preDepthCont1(:, recR), pind);
        [synth1, weight1] = ...
                nufftWarp(pyr1(:, :) .* AA1, ...
                depChg1(:, recR), -mem.depthClp1(:, recR), pind, filters, 2, 4, depPhsRatios, 1, minLev);
        pyr2 = occulsion(mem.pyr2(:, recR, c), depChg2(:, recR), mem.depth2(:, recR), preDepChg2(:, recR), mem.preDepthCont2(:, recR), pind);
        [synth2, weight2] = ...
                nufftWarp(pyr2(:, :) .* AA2, ...
                depChg2(:, recR), -mem.depthClp2(:, recR), pind, filters, 2, 4, depPhsRatios, -1, minLev);
        pyrRec = (synth1 .* weight1 + synth2 .* weight2) ./ (weight1 + weight2);
    else
        pyr2 = occulsion(mem.pyr2(:, recR, c), depChg2(:, recR), mem.depth2(:, recR), preDepChg2(:, recR), mem.preDepthCont2(:, recR), pind);
        pyrRec = nufftWarp(pyr2(:, :) .* AA2, ...
                depChg2(:, recR), -mem.depthClp2(:, recR), pind, filters, 2, 4, depPhsRatios, 1, minLev);
    end
end

function remapped = remapping(arr, method, origMinD, origMaxD, gama, fixZero)
    if method == 1,
        arr(arr < origMinD) = origMinD;
        arr(arr > origMaxD) = origMaxD;
        arr = (arr - origMinD) / (origMaxD - origMinD);
        arr = arr .^ gama;
        if ~fixZero,
            remapped = arr * origMaxD;
        else,
            zl = ((0 - origMinD) / (origMaxD - origMinD)) ^ gama;
            remapped = ((arr - zl) / (1-zl)) * origMaxD;
        end
    elseif method == 2,
        remapped = sign(arr) .* (abs(arr) .^ gama);
    else,
        remapped = arr;
    end
    remapped(isnan(remapped)) = 0;
end

function reduceArr = reduceToAvg(arr, siz, expand)
    reduceArr = zeros(siz);
    nR = size(arr, 1);
    nC = size(arr, 2);
    nSamp = ceil(expand * nR / siz(1));
    arrExp = zeros(nR+2*nSamp, nC);
    arrExp((1:nR) + nSamp, :) = arr;
    lSt = round(((1:siz(1)) - 0.5) * nR / siz(1) + 1 - nSamp/2);
    for dR = 0:nSamp-1,
        reduceArr = reduceArr + arrExp(lSt+dR+nSamp, :);
    end
    reduceArr = reduceArr / nSamp;
end

function reduceArr = reduce(fcn, arr, fac, ext)
    nR = size(arr, 1);
    reduceArr = arr(1:fac:nR, :);
    for ind = 2:fac,
        reduceArr = bsxfun(fcn, reduceArr, arr(ind:fac:nR, :));
    end
    for ind = 1-ext:0,
        reduceArr(2:nR/fac, :) = bsxfun(fcn, reduceArr(2:nR/fac, :), arr(ind+fac:fac:nR-fac, :));
    end
    for ind = fac+1:fac+ext,
        reduceArr(1:nR/fac-1, :) = bsxfun(fcn, reduceArr(1:nR/fac-1, :), arr(ind:fac:nR, :));
    end
end

function outFilters = get1DRadialFilters(dims, rVals, twidth)
    % Construct log_rad
    ctr = ceil((dims+0.5)/2);
    log_rad = abs(((1:dims)-ctr)./(dims/2));
    log_rad(ctr) = log_rad(ctr-1) / 4;

    N = max(size(rVals));

    [himask, lomaskPrev] = getRadialMaskPair(rVals(1),log_rad, twidth);
    mask = (1:dims)>ctr;
    outFilters{1} = himask.*mask;
    for k = 2:N
        %figure;
        [himask, lomask] = getRadialMaskPair(rVals(k),log_rad, twidth);
        outFilters{k} = himask.*lomaskPrev.*mask;
        lomaskPrev = lomaskPrev.*lomask;
        %plot(real(ifft(ifftshift(outFilters{k}))))
    end
    outFilters{N+1} = lomaskPrev;
end

function [himask, lomask] = getRadialMaskPair( r, log_rad, twidth)
    if false, %XXX
        n = length(log_rad);
        ctr = ceil((n + 0.5)/2);
        sig = r * (n / 2);
        rang = max(1, ceil(ctr - 2*sig)):min(n, floor(ctr + 2*sig));
        lomask = zeros(size(log_rad));
        lomask(rang) = exp(- (rang - ctr).^2 / (2*sig^2));
        himask = sqrt(1 - lomask.^2);
        return;
    end
    log_rad  = log2(log_rad)-log2(r);

    himask = log_rad;
    himask = clip(himask, -twidth, 0 );
    himask = himask * pi/(2*twidth);
    himask = abs(cos(himask));
    lomask = sqrt(1-himask.^2);
end

function ratios = getDepthPhaseRatio(filters)
    nLev = max(size(filters));
    nC = length(filters{1});
    for k = 1:nLev,
        rotated = ifftshift(filters{k});
        ratios{k} = nC / (2 * pi * imag(sum(rotated .* (0:nC-1) .* 1i) / sum(rotated)));
    end
end

function pind = buildBandIndices(filters, expand)
    nLev = max(size(filters));
    pind = [];
    sumIdx = 1;
    for k = 1:nLev
        curFilter = filters{k};
        indices = getIDXFromFilter1D(curFilter, expand);
        pind = [pind; sumIdx, numel(indices)];
        sumIdx = sumIdx + numel(indices);
    end
end

function indices =  pyrBandIndices(pind,band)
    indices = pind(band, 1):pind(band, 1) + pind(band, 2) - 1;
end

function pyr = buildSCFpyrGen1D(im, filters, expand)
    % Return pyramid in the usual format of a stack of column vectors
    imdft = fftshift(fft(im)); %DFT of image
    nLev = max(size(filters));
    n = length(im);
    pyr = [];
    for k = 1:nLev
        curFilter = filters{k};
        tempDFT = curFilter.*imdft; % Transform domain

        indices = getIDXFromFilter1D(curFilter, expand);
        nUnex = length(getIDXFromFilter1D(curFilter, 1));
        % shift to avoid cross-boundary wavelet
        tempDFT = tempDFT .* exp(2*pi*1i*((0:n-1) - floor(n/2))/nUnex);
        tempDFT = tempDFT(indices);
        curResult = ifft(ifftshift(tempDFT)) * length(indices) / nUnex;
        curResult = circshift(curResult, [0, floor(length(indices)/(2*nUnex))]);
        pyr = [pyr; curResult(:)];
    end
end

function pyr = buildSCFpyrGen1DLowpass(im, filters, expand, lev)
    % Return pyramid in the usual format of a stack of column vectors
    imdft = fftshift(fft(im)); %DFT of image
    nLev = max(size(filters));
    n = length(im);

    curFilter = filters{lev};
    [indices, lenO] = getIDXFromFilter1D(curFilter, expand);
    nUnex = length(getIDXFromFilter1D(curFilter, 1));

    curFilter = fftshift(fft(abs(ifft(ifftshift(curFilter)))));

    tempDFT = curFilter.*imdft; % Transform domain

    % shift to avoid cross-boundary wavelet
    tempDFT = tempDFT .* exp(2*pi*1i*((0:n-1) - floor(n/2))/nUnex);
    tempDFT = tempDFT(indices);
    curResult = ifft(ifftshift(tempDFT)) * length(indices) / nUnex;
    curResult = circshift(curResult, [0, floor(length(indices)/(2*nUnex))]);
    pyr = real(curResult(:));
end

function pyrNew = occulsion(pyr, mov, depth, oMov, oDepth, pind)
    pyrNew = pyr;
    nR = size(pyr, 2);
    nLev = size(pind, 1);
    n = pind(1, 2);
    for k = nLev-1:-1:1,
        indices = pyrBandIndices(pind, k);
        nL = length(indices);
        mt = 2*(mov(indices, :)*nL/n + repmat((0:nL-1)', [1, nR]));
        mtFloor = round(mt);
        oMt = 2*(oMov(indices, :)*nL/n + repmat((0:nL-1)', [1, nR]));
        oMtFloor = round(oMt);

        renderList = [reshape(repmat(1:nR, [nL, 1]), [nL*nR, 1]), ...
                      mtFloor(:), reshape(depth(indices, :), [nL*nR, 1]), ...
                      mt(:), reshape(repmat(1:nL, [nR, 1])', [nL*nR, 1]), ones(nL*nR, 1)];
        %renderList = [renderList; reshape(repmat(1:nR, [nL, 1]), [nL*nR, 1]), ...
        %              oMtFloor(:), reshape(oDepth(indices, :) + 4, [nL*nR, 1]), ...
        %              oMt(:), reshape(repmat(1:nL, [nR, 1])', [nL*nR, 1]), 2*ones(nL*nR, 1)];

        renderList = sortrows(renderList);
        renderList = renderList(size(renderList, 1):-1:1, :);
        [~, ia, ~] = unique(renderList(:, 1:2), 'rows', 'legacy');
        renderList = renderList(ia, :);

        nRend = size(renderList, 1);
        locs = sub2ind([nL, nR], renderList(:, 5), renderList(:, 1));
        disL = [(renderList(2:nRend, 4) - renderList(1:nRend-1, 4))/2; 1];
        disL(find(renderList(1:nRend-1, 3) <= renderList(2:nRend, 3))) = 1;
        disL(disL < 0) = 1;
        disL(disL > 1) = 1;
        disL2 = [(renderList(3:nRend, 4) - renderList(1:nRend-2, 4))/2; 1; 1];
        disL2(find(renderList(1:nRend-2, 3) <= renderList(3:nRend, 3))) = 1;
        disL2(disL2 < 0) = 1;
        disL2(disL2 > 1) = 1;
        disL = min(disL, disL2);
        disR = [1; (renderList(2:nRend, 4) - renderList(1:nRend-1, 4))/2];
        disR(find(renderList(1:nRend-1, 3) >= renderList(2:nRend, 3))+1) = 1;
        disR(disR < 0) = 1;
        disR(disR > 1) = 1;
        disR2 = [1; 1; (renderList(3:nRend, 4) - renderList(1:nRend-2, 4))/2];
        disR2(find(renderList(1:nRend-2, 3) >= renderList(3:nRend, 3))+2) = 1;
        disR2(disR2 < 0) = 1;
        disR2(disR2 > 1) = 1;
        disR = min(disR, disR2);
        dis = 2*(disL + disR) - 3;
        dis = max(dis, 0);
        w = 3 * dis .^ 2 - 2 * dis .^ 3;

        weight = zeros(nL, nR);
        mask = renderList(:, 6) == 1;
        weight(locs(mask)) = w(mask);
        %imshow(imresize(weight', [nR, n]));
        %figure;

        pyrNew(indices, :) = pyrNew(indices, :) .* weight;
    end
    %pause

    if 0
        if direction < 0,
            cum = zeros(nL, rSize);
            curMin = 1e9*ones(1, rSize);
            for iL = nL:-1:1,
                cum(iL, :) = curMin;
                curMin = min(curMin, mt(iL, :));
            end
            %fac = exp(-2*(max((mt - cum)/m + 1, 0)));
            fac = (mt - cum)/m + 1;
            fac = 1 - min(max(fac * 4, 0), 1);
            fac = -2 * fac .^ 3 + 3 * fac .^ 2;
            pyr(indices, rInd) = pyr(indices, rInd) .* fac;
        else,
            cum = zeros(nL, rSize);
            curMax = -1e9*ones(1, rSize);
            for iL = 1:nL,
                cum(iL, :) = curMax;
                curMax = max(curMax, mt(iL, :));
            end
            %fac = exp(-2*(max((cum - mt)/m + 1, 0)));
            fac = (cum - mt)/m + 1;
            fac = 1 - min(max(fac * 4, 0), 1);
            fac = -2 * fac .^ 3 + 3 * fac .^ 2;
            pyr(indices, rInd) = pyr(indices, rInd) .* fac;
        end
    end
end

function [resPyr, weight] = nufftWarp(pyr, mov, depth, pind, filters, m, Q, depPhsRatios, direction, minLev)
    %http://www.cscamm.umd.edu/programs/fam04/qing_liu_fam04.pdf
    nLev = max(size(filters));
    nR = size(pyr, 2);
    n = length(filters{1});
    resPyr = zeros(size(pyr));
    weight = ones(size(pyr));
    for rSt = 1:50:nR,
        rEd = min(rSt+49, nR);
        rSize = rEd-rSt+1;
        rInd = rSt:rEd;
        for k = nLev:-1:minLev,
            curFilter = filters{k};
            indices = pyrBandIndices(pind, k);

            %if k == nLev,
            %    mov(indices, :) = 0;
            %end
            if false%k == nLev,
                %resPyr(indices, rInd) = pyr(indices, rInd);
                nL = pind(k, 2);
                for rId = rSt:rEd,
                    rLoc = mov(indices, rId)*nL/n + (0:nL-1)';
                    resPyr(indices, rId) = interp1([rLoc-nL; rLoc; rLoc+nL], repmat(pyr(indices, rId), [3, 1]), (0:nL-1)', 'pchip');
                end
                continue;
            end
            nL = pind(k, 2);
            q = floor(min(Q, nL-1)/2)*2;

            %q, nL
            F = zeros(q+1);
            for iQ = 1:q+1,
                for jQ = 1:q+1,
                    if mod(iQ - jQ, m*nL) == 0,
                        F(iQ, jQ) = nL;
                    else,
                        zINo2 = exp(1i*pi*(jQ-iQ)/m);
                        zI = exp(1i*2*pi*(jQ-iQ)/(m*nL));
                        F(iQ, jQ) = (1/zINo2 - zINo2) / (1 - zI);
                    end
                end
            end
            F = inv(F);
            A = zeros(q+1, nL, rSize);
            mt = m*(mov(indices, rInd)*nL/n + repmat((0:nL-1)', [1, rSize]));
            mtFloor = round(mt);
            stack = zeros(nL, 1);
            W = ones(nL, rSize);

            if 0
                for rId = rInd,
                    if direction < 0,
                        mask = mt(:, rId-rSt+1) < cum(:, rId-rSt+1);
                    else,
                        mask = cum(:, rId-rSt+1) < mt(:, rId-rSt+1);
                    end
                    locs = find(mask);

                    nLoc = size(locs, 1);
                    dis = mov(indices, rId);
                    dis = dis(locs);

                    movs = mt(mask, rId-rSt+1) / m;
                    movEx = zeros(nLoc+2, 1);
                    movEx(2:nLoc+1) = movs;
                    movEx(1) = movs(1);
                    movEx(nLoc+2) = movs(nLoc);
                    dif = movEx(3:nLoc+2) - movEx(1:nLoc);
                    dif(2:nLoc-1) = dif(2:nLoc-1)/2;

                    w1 = (abs(dis) + 1e-3);
                    w2 = (cos((pi / 4) * min(2, max(0, dif - 1))));
                    wr = interp1(movs, [w1, w2], (0:nL-1)', 'pchip');
                    weight(indices, rId) = wr(:, 2) ./ max(wr(:, 1), 1e-4);
                end
            end

            for iQ = -q/2:q/2,
                frac = mt - mtFloor - iQ;
                tmp1 = -1i*(sin((pi/m)*(frac-0.5))./(1-exp(1i*2*pi*(frac-0.5)/(nL*m))));
                tmp2 = -1i*(sin((pi/m)*(frac+0.5))./(1-exp(1i*2*pi*(frac+0.5)/(nL*m))));
                tmp1(abs(frac-0.5)<1e-10) = nL/2;
                tmp2(abs(frac+0.5)<1e-10) = nL/2;

                A(iQ+q/2+1, :, :) = tmp1+tmp2;
            end

            mtFloor = mod(mtFloor, nL*m);

            X = reshape(F * reshape(A, [q+1, nL*rSize]), [q+1, nL, rSize]);
            sig = zeros(m*nL+q, rSize);
            sig2 = zeros(m*nL+q, rSize);
            for iQ = 1:q+1,
                subs = [reshape(mtFloor+iQ, [nL*rSize, 1]), reshape(repmat(1:rSize, [nL, 1]), [nL*rSize, 1])];
                sig = sig + accumarray(subs, reshape(pyr(indices, rInd) .* W .* reshape(X(iQ, :, :), [nL, rSize]),...
                        [nL*rSize, 1]), [m*nL+q, rSize]);
            end
            if q/2 <= m * nL,
                sig(1+m*nL:q/2+m*nL, :) = sig(1+m*nL:q/2+m*nL, :) + sig(1:q/2, :);
                sig(q/2+1:q, :) = sig(q/2+1:q, :) + sig(q/2+1+m*nL:q+m*nL, :);
            else,
                for iQ = 1:q/2,
                    sig(mod(iQ + m*L - mod(q/2, m*L) - 1, m*L) + q/2 + 1, :) = ...
                            sig(mod(iQ + m*L - mod(q/2, m*L) - 1, m*L) + q/2 + 1, :) + sig(iQ, :);
                    sig(mod(iQ - 1, m*L) + q/2 + 1, :) = sig(mod(iQ - 1, m*L) + q/2 + 1, :) + sig(q/2+iQ+m*nL, :);
                end
            end
            sig = sig(q/2+1:q/2+m*nL, :);
            for iR = 1:rSize,
                fftSig = fftshift(fft(sig(:, iR)));
                fftSig = fftSig ./ cos(pi*((0:m*nL-1)-m*nL/2)/(m*nL))';
                resPyr(indices, rSt+iR-1) = ifft(ifftshift(fftSig(m*nL/2-floor(nL/2)+(1:nL))));
            end
        end
    end
end

function res = reconSCFpyrGen1D(pyr, pind, filters, minLev, maxLev)
    nLev = max(size(filters));
    n = length(filters{1});
    imdft = zeros(1, n); %DFT of image
    for k = minLev:maxLev
        curFilter = filters{k};
        tempDFT = fftshift(fft(2*real(pyr(pyrBandIndices(pind,k))'))); % Transform domain
        if k < nLev,
            curFilter = curFilter + [0, rot90(curFilter(2:n),2)];
        else
            tempDFT = tempDFT / 2;
        end
        indices = getIDXFromFilter1D(curFilter, 1);
        imdft(indices) = imdft(indices) + tempDFT.*curFilter(indices)...
                .*exp(-2*pi*1i*(indices-1-floor(n/2))/length(indices));
    end
    res =  real(ifft(ifftshift(imdft)));
end

function [ out, lenO ] = getIDXFromFilter1D(filt, expand)
    aboveZero = filt>1e-10;
    aboveZero = logical(aboveZero + fliplr(aboveZero));
    n = length(filt);

    idx1 = 1:length(filt);
    idx1 = idx1(aboveZero);
    lenO = max(idx1) - min(idx1) + 1;
    len = lenO;

    fac = 1;
    while fac * 2 <= expand && len * fac * 2 <= n,
        fac = fac * 2;
    end
    len = len * fac;
    ctr = ceil((n + 1) / 2);
    if mod(len, 2) == 0,
        out = ctr - len/2:ctr + len/2-1;
    else,
        out = ctr - (len-1)/2:ctr + (len-1)/2;
    end
end

function filtered = bilateral(img, guid, wgt, sig_h, sig_v, sig_c)
    nR = size(img, 1);
    nC = size(img, 2);
    nCh = size(img, 3);
    nCg = size(guid, 3);

    hEx = ceil(sig_h * 2);
    vEx = ceil(sig_v * 2);

    imgEx = zeros(nR + vEx*2, nC + hEx*2, nCh);
    imgEx(vEx + (1:nR), hEx + (1:nC), :) = img;
    guidEx = zeros(nR + vEx*2, nC + hEx*2, nCg);
    guidEx(vEx + (1:nR), hEx + (1:nC), :) = guid;
    frameEx = zeros(nR + vEx*2, nC + hEx*2);
    frameEx(vEx + (1:nR), hEx + (1:nC)) = wgt;

    filtered = zeros(nR, nC, nCh);
    weight = zeros(nR, nC, nCh);

    for dR = -vEx:vEx,
        for dC = -hEx:hEx,
            dW = exp(-0.5*((dR^2)/(sig_v^2) + (dC^2)/(sig_h^2) + ...
                    sum((guidEx(vEx+dR+(1:nR), hEx+dC+(1:nC), :) - guid(:, :, :)).^2, 3) / (sig_c^2)));
            dW = dW .* frameEx(vEx+dR+(1:nR), hEx+dC+(1:nC));
            dW = repmat(dW, [1, 1, nCh]);
            filtered = filtered + imgEx(vEx+dR+(1:nR), hEx+dC+(1:nC), :) .* dW;
            weight = weight + dW;
        end
    end

    filtered = filtered ./ weight;
end
function result = getShiftWeight(dR, dC, cSt, cEd, nR, sigR, sigC)
    result = exp(-0.5 * (dR*dR/(sigR*sigR)+dC*dC/(sigC*sigC))) * ones(nR, cEd-cSt+1);
end
function result = getShiftWeightConst(dR, dC, cSt, cEd, nR)
    result = ones(nR, cEd-cSt+1);
end
