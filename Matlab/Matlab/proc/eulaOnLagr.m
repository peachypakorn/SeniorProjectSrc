function [synth, depthLo, depthRo, phsCmpO, vis] = eulaOnLagr(img1, img2, depthLo, depthRo, testParam, runParam, noEula, phsCmpO)
    rho = runParam.rho;
    ht = runParam.ht;
    nA = runParam.nA;
    dA = runParam.dA;
    dsrLagr = runParam.eurLagrResize;
    maxDLagr = ceil(testParam.maxDepth * dsrLagr);
    vis = struct();
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
    if isfield(runParam, 'expPrealign'),
        expPrealign = runParam.expPrealign;
    else,
        expPrealign = false;
    end

    nR = size(img1, 1);
    nC = size(img1, 2);
    nCol = size(img1, 3);
    maxAmp = (nA-1)*dA;
    sigma_v_fac = 2;

    if false,
        % XXX force to power of two
        nCf = round(2^round(log(nC) / log(2)));
        img1 = imresize(img1, [nR, nCf]);
        img2 = imresize(img2, [nR, nCf]);
        tmp = nCf;
        nCf = nC;
        nC = tmp;
        synthT = zeros([ nR, nC, nCol,  nA*2]);
    end

    if isscalar(depthLo) || isscalar(depthRo),
        [lagrDepthL, lagrDepthR] = lagrDepthEstimate(imresize(img1, dsrLagr), imresize(img2, dsrLagr), maxDLagr, dReg);
        depthLo = lagrDepthL;
        depthRo = lagrDepthR;
        dsrLagr = size(lagrDepthL, 2) / size(img1, 2);
    else
        lagrDepthL = depthLo;
        lagrDepthR = depthRo;
        dsrLagr = size(lagrDepthL, 2) / size(img1, 2);
    end

    sizFilt = 2.^[0:-1:-ht];
    filters = get1DRadialFilters(nC, sizFilt, 1);
    nLev = max(size(filters));
    depPhsRatios = getDepthPhaseRatio(filters);

    pind = buildBandIndices(filters, 1);
    nF = pind(nLev, 1) + pind(nLev, 2) - 1;

    pindRef = buildBandIndices(filters, 16);
    nFRef = pindRef(nLev, 1) + pindRef(nLev, 2) - 1;

    pyr1 = zeros([nF nR nCol]);
    pyr2 = zeros([nF nR nCol]);

    pyr1Ref = zeros([nFRef nR nCol]);
    pyr2Ref = zeros([nFRef nR nCol]);
    refInd1 = ones([nF nR nCol]);
    refInd2 = ones([nF nR nCol]);
    preDepth1 = zeros([nF nR]);
    preDepth2 = zeros([nF nR]);
    preDepthCont1 = zeros([nF nR]);
    preDepthCont2 = zeros([nF nR]);

    for iR=1:nR,
        for iCol = 1:nCol,
            pyr1(:, iR, iCol) = buildSCFpyrGen1D(img1(iR, :, iCol), filters, 1);
            pyr1Ref(:, iR, iCol) = buildSCFpyrGen1D(img1(iR, :, iCol), filters, 16);
            pyr2(:, iR, iCol) = buildSCFpyrGen1D(img2(iR, :, iCol), filters, 1);
            pyr2Ref(:, iR, iCol) = buildSCFpyrGen1D(img2(iR, :, iCol), filters, 16);
        end
    end;

    vis.pyr1 = cell(nLev, 1);
    vis.pyr2 = cell(nLev, 1);
    vis.layer1 = cell(nLev, 1);
    vis.layer2 = cell(nLev, 1);
    vis.accum1 = cell(nLev, 1);
    vis.accum2 = cell(nLev, 1);
    for iLev = nLev:-1:1,
        ind = pyrBandIndices(pind, iLev);
        vis.pyr1{iLev} = pyr1(ind, 1, 1);
        vis.pyr2{iLev} = pyr2(ind, 1, 1);

        vis.layer1{iLev} = reconSCFpyrGen1D(pyr1(:, 1, 1), pind, filters, iLev, iLev);
        vis.accum1{iLev} = reconSCFpyrGen1D(pyr1(:, 1, 1), pind, filters, iLev, nLev);

        vis.layer2{iLev} = reconSCFpyrGen1D(pyr2(:, 1, 1), pind, filters, iLev, iLev);
        vis.accum2{iLev} = reconSCFpyrGen1D(pyr2(:, 1, 1), pind, filters, iLev, nLev);
    end

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

        if ~expPrealign,
            levDepthL = imresize(lagrDepthL, [nR, pind(iLev, 2)*8])';
            sigma_v = sigma_v_fac;
            h = fspecial('gaussian', [1, ceil(sigma_v*2)*2+1], sigma_v+1e-10);
            levDepthL = imfilter(levDepthL, h, 'symmetric');
            levDepthL = reduce(@min, levDepthL, 8, 8);
            preDepthCont1(ind, :) = levDepthL/dsrLagr;

            levDepthR = imresize(lagrDepthR, [nR, pind(iLev, 2)*8])';
            levDepthR = imfilter(levDepthR, h, 'symmetric');
            levDepthR = reduce(@min, levDepthR, 8, 8);
            preDepthCont2(ind, :) = levDepthR/dsrLagr;
        else,
            levDepthL = imresize(lagrDepthL, [nR, nC]);
            levDepthR = imresize(lagrDepthR, [nR, nC]);
            sigma_v = sigma_v_fac;
            h = fspecial('gaussian', [ceil(sigma_v*2)*2+1, 1], sigma_v+1e-10);
            levDepthLO = imfilter(levDepthL, h, 'symmetric');
            levDepthRO = imfilter(levDepthR, h, 'symmetric');
            levDepthL = zeros([pind(iLev, 2), nR]);
            levDepthR = zeros([pind(iLev, 2), nR]);

            imgGuide = zeros(nR, nC);
            imgGuidePos = zeros(nR, nC);
            imgGuideNeg = zeros(nR, nC);

            if false,
                for iR = 1:nR,
                    if nCol == 3,
                        bw = pyr1(:, iR, 1) * 0.3 + pyr1(:, iR, 2) * 0.6 + pyr1(:, iR, 3) * 0.1;
                    else,
                        bw = pyr1(:, iR, 1);
                    end
                    imgGuide(iR, :) = reconSCFpyrGen1D(bw, pind, filters, iLev, iLev) ;
                end
            else,
                if nCol == 3,
                    imBW = img1(:, :, 1) * 0.3 + img1(:, :, 2) * 0.6 + img1(:, :, 3) * 0.1;
                else,
                    imBW = img1;
                end
                sigma_bp = 2^(iLev-1);
                h1 = fspecial('gaussian', [1, ceil(sigma_bp*2)*2+1], sigma_bp);
                h2 = fspecial('gaussian', [1, ceil(sigma_bp*4)*2+1], sigma_bp*2);
                imgGuide = imfilter(imBW, h1, 'symmetric') - imfilter(imBW, h2, 'symmetric');
            end


            imgGuidePos(imgGuide > 0) = imgGuide(imgGuide > 0) + 1e-5;
            imgGuideNeg(imgGuide < 0) = -imgGuide(imgGuide < 0) + 1e-5;
            for iR = 1:nR,
                lowpassPos = buildSCFpyrGen1DLowpass(imgGuidePos(iR, :) .* levDepthLO(iR, :), filters, 1, iLev);
                lowpassPosW = buildSCFpyrGen1DLowpass(imgGuidePos(iR, :), filters, 1, iLev);
                lowpassNeg = buildSCFpyrGen1DLowpass(imgGuideNeg(iR, :) .* levDepthLO(iR, :), filters, 1, iLev);
                lowpassNegW = buildSCFpyrGen1DLowpass(imgGuideNeg(iR, :), filters, 1, iLev);
                reduced = bsxfun(@min, lowpassPos ./ lowpassPosW, lowpassNeg ./ lowpassNegW);
                preDepthCont1(ind, iR) = preDepthCont1(ind, iR) + reduced / (dsrLagr);
                levDepthL(:, iR) = levDepthL(:, iR) + reduced;
            end
            out = zeros([nR, nC, 3]);
            out(:, :, 1) = imgGuidePos;
            out(:, :, 2) = imgGuideNeg;
            out(:, :, 3) = imgGuideNeg;
            %imwrite(imresize(repmat((levDepthL'*4 + 128) / 255, [1, 1, 3]), [nR, nC]), sprintf('tmp/d%03d.jpg', iLev));
            %imwrite(out * 4, sprintf('tmp/%03d.jpg', iLev));

            imgGuide = zeros(nR, nC);
            imgGuidePos = zeros(nR, nC);
            imgGuideNeg = zeros(nR, nC);

            if false,
                for iR = 1:nR,
                    if nCol == 3,
                        bw = pyr2(:, iR, 1) * 0.3 + pyr2(:, iR, 2) * 0.6 + pyr2(:, iR, 3) * 0.1;
                    else,
                        bw = pyr2(:, iR, 1);
                    end
                    imgGuide(iR, :) = reconSCFpyrGen1D(bw, pind, filters, iLev, iLev) ;
                end
            else
                if nCol == 3,
                    imBW = img2(:, :, 1) * 0.3 + img2(:, :, 2) * 0.6 + img2(:, :, 3) * 0.1;
                else,
                    imBW = img2;
                end
                sigma_bp = 2^(iLev-2);
                h1 = fspecial('gaussian', [1, ceil(sigma_bp*2)*2+1], sigma_bp);
                h2 = fspecial('gaussian', [1, ceil(sigma_bp*4)*2+1], sigma_bp*2);
                imgGuide = imfilter(imBW, h1, 'symmetric') - imfilter(imBW, h2, 'symmetric');
            end

            imgGuidePos(imgGuide > 0) = imgGuide(imgGuide > 0) + 1e-5;
            imgGuideNeg(imgGuide < 0) = -imgGuide(imgGuide < 0) + 1e-5;

            for iR = 1:nR,
                lowpassPos = buildSCFpyrGen1DLowpass(imgGuidePos(iR, :) .* levDepthRO(iR, :), filters, 1, iLev);
                lowpassPosW = buildSCFpyrGen1DLowpass(imgGuidePos(iR, :), filters, 1, iLev);
                lowpassNeg = buildSCFpyrGen1DLowpass(imgGuideNeg(iR, :) .* levDepthRO(iR, :), filters, 1, iLev);
                lowpassNegW = buildSCFpyrGen1DLowpass(imgGuideNeg(iR, :), filters, 1, iLev);
                reduced = bsxfun(@min, lowpassPos ./ lowpassPosW, lowpassNeg ./ lowpassNegW);
                preDepthCont2(ind, iR) = preDepthCont2(ind, iR) + reduced / (dsrLagr);
                levDepthR(:, iR) = levDepthR(:, iR) + reduced;
            end
        end

        preDepth1(ind, :) = repmat((-0.5 + (1:nO)')*(nC/nO), [1, nR]) - levDepthL/dsrLagr;
        preDepth1(ind, :) = clip(round(preDepth1(ind, :) * (nRef / nC) - refStart), 0, nRef-1);
        refInd1(ind, :, :) = repmat(preDepth1(ind, :) + pindRef(iLev, 1), [1, 1, nCol])+ ...
                             repmat((0:nR-1)*nFRef, [length(ind), 1, nCol]) + ...
                             repmat(reshape((0:nCol-1)*nFRef*nR, [1, 1, nCol]), [length(ind), nR, 1]);
        preDepth1(ind, :) = (preDepth1(ind, :) + refStart) * (nC / nRef) - ...
                repmat((-0.5 + (1:nO)')*(nC/nO), [1, nR]);

        preDepth2(ind, :) = repmat((-0.5 + (1:nO)')*(nC/nO), [1, nR]) + levDepthR/dsrLagr;
        preDepth2(ind, :) = clip(round(preDepth2(ind, :) * (nRef / nC) - refStart), 0, nRef-1);
        refInd2(ind, :, :) = repmat(preDepth2(ind, :) + pindRef(iLev, 1), [1, 1, nCol]) + ...
                             repmat((0:nR-1)*nFRef, [length(ind), 1, nCol]) + ...
                             repmat(reshape((0:nCol-1)*nFRef*nR, [1, 1, nCol]), [length(ind), nR, 1]);
        preDepth2(ind, :, :) = repmat((-0.5 + (1:nO)')*(nC/nO), [1, nR]) - ...
                (preDepth2(ind, :) + refStart) * (nC / nRef);
    end

    %figure
    %hold on;
    %val = pyr1(pyrBandIndices(pind, 5), floor(nR/2), 2);
    %plot(angle(val(2:pind(5, 2)) ./ val(1:pind(5, 2)-1)), 'red');
    %plot(abs(val));
    %figure
    %hold on;
    %val = pyr1(pyrBandIndices(pind, 4), floor(nR/2), 2);
    %plot(angle(val(2:pind(4, 2)) ./ val(1:pind(4, 2)-1)), 'red');
    %plot(abs(val));
    %figure

    %% Computing Filter
    d1 = angle(pyr1);
    d2 = angle(pyr2);
    d1Ref = angle(pyr1Ref);
    d2Ref = angle(pyr2Ref);
    phsDif1 = angle(pyr2Ref(refInd1)) - angle(pyr1);
    phsDif2 = angle(pyr2) - angle(pyr1Ref(refInd2));
    weight1 = bsxfun(@min, abs(pyr2Ref(refInd1)), abs(pyr1));
    weight2 = bsxfun(@min, abs(pyr2), abs(pyr1Ref(refInd2)));

    sigma_v = sigma_v_fac;
    ind = pyrBandIndices(pind, nLev);
    h = fspecial('gaussian', [ceil(sigma_v*2)*2+1, 1], sigma_v+1e-10);
    tmpDepthL = -imfilter(imresize(lagrDepthL, [nR, pind(nLev, 2)*4]), h, 'symmetric')'/dsrLagr;
    preDepth1(ind, :) = reduce(@max, tmpDepthL, 4, 2);
    tmpDepthR = -imfilter(imresize(lagrDepthR, [nR, pind(nLev, 2)*4]), h, 'symmetric')'/dsrLagr;
    preDepth2(ind, :) = reduce(@max, tmpDepthR, 4, 2);
    phsDif1(ind, :, :) = 0;
    phsDif2(ind, :, :) = 0;
    weight1(ind, :, :) = 1;
    weight2(ind, :, :) = 1;

    for k = nLev-1:-1:1,
        ind = pyrBandIndices(pind,k);
        phsDif1(ind, :, :) = repmat(-preDepth1(ind, :) / depPhsRatios{k}, [1, 1, nCol]) + phsDif1(ind, :, :);
        phsDif2(ind, :, :) = repmat(-preDepth2(ind, :) / depPhsRatios{k}, [1, 1, nCol]) + phsDif2(ind, :, :) ;
    end

    phsDif1 = exp(phsDif1*1i);
    phsDif2 = exp(phsDif2*1i);
    wPhs1 = ones([nF, nR])*1e-20;
    wPhs2 = ones([nF, nR])*1e-20;
    sumPhs1 = zeros([nF, nR]);
    sumPhs2 = zeros([nF, nR]);

    for k = nLev-1:-1:1,
        ind = pyrBandIndices(pind,k);
        cLen = (2^max(0,k-2));
        sigma_v = sigma_v_fac;
        for c = 1:nCol,
            for d = 0:sigma_v*2,
                wTmp = squeeze(weight1(ind, 1+d:nR, c)) * exp(-(d).^2 / (2 * sigma_v^2))';
                sumPhs1(ind, 1:nR-d) = sumPhs1(ind, 1:nR-d) + phsDif1(ind, 1+d:nR, c) .* wTmp;
                wPhs1(ind, 1:nR-d) = wPhs1(ind, 1:nR-d) + wTmp;

                wTmp = squeeze(weight2(ind, 1+d:nR, c)) * exp(-(d).^2 / (2 * sigma_v^2))';
                sumPhs2(ind, 1:nR-d) = sumPhs2(ind, 1:nR-d) + phsDif2(ind, 1+d:nR, c) .* wTmp;
                wPhs2(ind, 1:nR-d) = wPhs2(ind, 1:nR-d) + wTmp;

                if d == 0,
                    continue;
                end

                wTmp = squeeze(weight1(ind, 1:nR-d, c)) * exp(-(d).^2 / (2 * sigma_v^2))';
                sumPhs1(ind, 1+d:nR) = sumPhs1(ind, 1+d:nR) + phsDif1(ind, 1:nR-d, c) .* wTmp;
                wPhs1(ind, 1+d:nR) = wPhs1(ind, 1+d:nR) + wTmp;

                wTmp = squeeze(weight2(ind, 1:nR-d, c)) * exp(-(d).^2 / (2 * sigma_v^2))';
                sumPhs2(ind, 1+d:nR) = sumPhs2(ind, 1+d:nR) + phsDif2(ind, 1:nR-d, c) .* wTmp;
                wPhs2(ind, 1+d:nR) = wPhs2(ind, 1+d:nR) + wTmp;
            end;
        end
    end
    phsDif1 = angle(sumPhs1 ./ wPhs1);
    phsDif2 = angle(sumPhs2 ./ wPhs2);
    ind = pyrBandIndices(pind, nLev);
    phsDif1(ind, :) = 0;
    phsDif2(ind, :) = 0;

    for k = nLev-1:-1:1,
        ind = pyrBandIndices(pind,k);
        phsDif1(ind, :) = mod(phsDif1(ind, :) + preDepth1(ind, :) / depPhsRatios{k} + pi, 2*pi) - pi;
        phsDif2(ind, :) = mod(phsDif2(ind, :) + preDepth2(ind, :) / depPhsRatios{k} + pi, 2*pi) - pi;
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

    if maxDepth > 0,
        depthClp1 = maxDepth * atan(depth1 / maxDepth);
        depthClp2 = maxDepth * atan(depth2 / maxDepth);
        %depthClp = clip(depth, -maxDepth, maxDepth);
    else
        depthClp1 = depth1;
        depthClp2 = depth2;
    end

    phsDifAA1 = zeros([nF, nR]);
    phsDifAA2 = zeros([nF, nR]);
    for k = nLev-1:-1:1,
        ind = pyrBandIndices(pind,k);
        phsDifAA1(ind, :) = dA * depthClp1(ind, :) / depPhsRatios{k};
        phsDifAA2(ind, :) = dA * depthClp2(ind, :) / depPhsRatios{k};
    end

    if recR == 0,
        recR = 1:nR;
        synth = zeros([ nR, nC, nCol,  nA*2]);
    else
        synth = zeros([ length(recR) nC, nCol,  nA*2]);
    end
    vis.mov = cell(2*nA, 1);
    mulCorr = ones([nF length(recR) nA*2+appWrap*2]);
    for c = 1:nCol,
        indices = pyrBandIndices(pind, nLev);
        %tmp = (pyr1(indices, recR, c) + pyr2(indices, recR, c)) / 2;
        %pyr1(indices, recR, c)  = tmp;
        %pyr2(indices, recR, c)  = tmp;
        %% Reconstructing
        pyrRec = zeros(nF, length(recR), 2*nA+2*appWrap);
        if expShear,
            aRange = [nA + 0.5, 1:nA+appWrap];
        else,
            aRange = 1:nA+appWrap;
        end
        for iA = aRange,
            amp = -0.5 + (iA - 0.5)*dA;

            depChg11 = -depthClp1 * (0.5+amp) + depth1/2;
            depChg12 = -depthClp1 * (0.5+amp) - depth1/2;
            depChg21 = -depthClp2 * (0.5+amp) + depth2/2;
            depChg22 = -depthClp2 * (0.5+amp) - depth2/2;
            phsChg11 = zeros(size(depChg11));
            phsChg12 = zeros(size(depChg12));
            phsChg21 = zeros(size(depChg21));
            phsChg22 = zeros(size(depChg22));
            fcc11 = abs(depChg11);
            fcc12 = abs(depChg12);
            fcc21 = abs(depChg21);
            fcc22 = abs(depChg22);
            ind1 = fcc11 + fcc22 == 0;
            ind2 = fcc21 + fcc12 == 0;
            fcc11(ind1) = abs(amp*1e-10);
            fcc21(ind2) = abs(amp*1e-10);
            fcc12(ind2) = abs((1+amp)*1e-10);
            fcc22(ind1) = abs((1+amp)*1e-10);
            mixFactor1 = fcc22 ./ (fcc11 + fcc22);
            mixFactor2 = fcc12 ./ (fcc12 + fcc21);
            %mixFactor1 = ones(size(fcc11));
            %mixFactor2 = ones(size(fcc11));
            for k = nLev:-1:1,
                ind = pyrBandIndices(pind,k);
                phsChg11(ind, :) = depChg11(ind, :) / depPhsRatios{k};
                phsChg12(ind, :) = depChg12(ind, :) / depPhsRatios{k};
                phsChg21(ind, :) = depChg21(ind, :) / depPhsRatios{k};
                phsChg22(ind, :) = depChg22(ind, :) / depPhsRatios{k};
            end

            if iA > nA - appWrap && iA ~= nA + 0.5,
                AA1 = exp( -( phsDifAA1(:, recR) * rho).^2 / 2);
                AA2 = exp( -( phsDifAA2(:, recR) * rho).^2 / 2);
                sig = sigWrap * ones(nF, length(recR));
                for k = nLev:-1:1,
                    ind = pyrBandIndices(pind,k);
                    sig(ind, :) = sig(ind, :) / (2^(nLev - k));%, appWrap/2);
                end
                sig(sig > appWrap * dA/2) = appWrap * dA/2;
                %pref = repmat(0.5 + 0.5 * erf((nA + 0.5 - iA) * dA./ (sqrt(2) * sig(:, 1))), [1, length(recR)]);
                AA1 = 0.5 * AA1 .* (1 + erfz((-dA * (iA - nA - 0.5)  - phsDifAA1(:, recR) * rho^2 * dA * i) ...
                        ./ sqrt(2 * (rho^2 * dA^2) + sig(:, :).^2)));
                %AA1 = AA1 ./ pref;
                AA2 = 0.5 * AA2 .* (1 + erfz((-dA * (iA - nA - 0.5)  + phsDifAA2(:, recR) * rho^2 * dA * i) ...
                        ./ sqrt(2 * (rho^2 * dA^2) + sig(:, :).^2)));
                %AA2 = AA2 ./ pref;

            elseif expSinc,
                AA1 = phsDifAA1(:, recR) * rho < pi/2;
                AA2 = phsDifAA2(:, recR) * rho < pi/2;
            else,
                AA1 = exp( -( phsDifAA1(:, recR) * rho).^2 / 2);
                AA2 = exp( -( phsDifAA2(:, recR) * rho).^2 / 2);
            end
            if iA == nA + 0.5,
                if c > 1,
                    continue;
                end
                pyrCmp1 = nufftWarp((pyr1(:, recR, 1) * 0.3 + pyr1(:, recR, 2) * 0.6 + pyr1(:, recR, 3) * 0.1) .* AA1, ...
                                -depChg11(:, recR), -depthClp1(:, recR), pind, filters, 2, 8, depPhsRatios);

                pyrCmp2 = nufftWarp((pyr2(:, recR, 1) * 0.3 + pyr2(:, recR, 2) * 0.6 + pyr2(:, recR, 3) * 0.1) .* AA2, ...
                                depChg21(:, recR), -depthClp2(:, recR), pind, filters, 2, 8, depPhsRatios);

                keyLev = runParam.shearLev;
                phsCmpA = 0;
                for curLev = keyLev:-1:runParam.shearMinLev,
                    po2 = 2^(keyLev - curLev);
                    gR = min(runParam.shearGrid(1) * po2, length(recR));
                    gC = runParam.shearGrid(2) * po2;
                    inds = pyrBandIndices(pind, curLev);
                    nL = length(inds);
                    nRR = length(recR);
                    phsCont = (nC / nL) / depPhsRatios{curLev};
                    imCmp1 = pyrCmp1(inds, :) .* repmat(exp(-i*(1:nL)'*phsCont), [1, nRR]);
                    imCmp2 = pyrCmp2(inds, :) .* repmat(exp(-i*(1:nL)'*phsCont), [1, nRR]);
                    imCmp1 = imresize(imCmp1, [gC, nRR]);
                    imCmp2 = imresize(imCmp2, [gC, nRR]);
                    dif = imCmp1 ./ imCmp2;
                    dif = dif .* bsxfun(@min, abs(imCmp1), abs(imCmp2)) ./ abs(dif);
                    phsCmp = angle(imresize(dif, [gC, gR])) / (2*nA * po2);
                    if ~isscalar(phsCmpA),
                        phsCmpA = imresize(phsCmpA, [gC, gR]);
                        phsCmpA = phsCmp + (pi / (nA*po2)) * round((phsCmpA - phsCmp) / (pi / (nA*po2)));
                    else
                        phsCmpA = phsCmp;
                    end
                end
                if ~isscalar(phsCmpO),
                        phsCmpA = phsCmpA + (pi / nA) * round((phsCmpO - phsCmpA) / (pi / nA));
                        phsCmpA = phsCmpO * 0.9 + phsCmpA * 0.1;
                end
                phsCmpO = phsCmpA;
                phsCmp = phsCmpA;
                for iA2 = 1:2*nA+2*appWrap,
                    for iLev = 1:nLev-1,
                        inds = pyrBandIndices(pind, iLev);
                        imCorr = exp(i * (depPhsRatios{keyLev} / depPhsRatios{iLev}) * phsCmp * (iA2-nA-appWrap-0.5));
                        mulCorr(inds, :, iA2) = imresize(imCorr, [length(inds), nRR]);
                    end
                end
                mulCorr = mulCorr ./ abs(mulCorr);

            else,
                if c == 1 && iA <= nA,
                    [tmp, movVis] = nufftWarp(pyr1(:, 1:1, c) .* AA1(:, 1:1), ...
                                -depChg11(:, 1:1), -depthClp1(:, 1:1), pind, filters, 2, 8, depPhsRatios);
                    vis.mov{nA+1-iA} = movVis;
                    [tmp, movVis] = nufftWarp(pyr2(:, 1:1, c) .* AA2(:, 1:1), ...
                                depChg21(:, 1:1), -depthClp2(:, 1:1), pind, filters, 2, 8, depPhsRatios);
                    vis.mov{nA+iA} = movVis;
                end

                pyrRec(:, :, nA+appWrap+1-iA) = ...
                        nufftWarp(pyr1(:, recR, c) .* AA1, ...
                                -depChg11(:, recR), -depthClp1(:, recR), pind, filters, 2, 8, depPhsRatios); %.* mixFactor1(:, recR) +...
                        %nufftWarp(pyr2(:, recR, c) .* AA2, ...
                        %        -depChg22(:, recR), -depthClp2(:, recR), pind, filters, 2, 8) .* (1-mixFactor1(:, recR));

                pyrRec(:, :, nA+appWrap+iA) = ...
                        nufftWarp(pyr2(:, recR, c) .* AA2, ...
                                depChg21(:, recR), -depthClp2(:, recR), pind, filters, 2, 8, depPhsRatios);% .* mixFactor2(:, recR) +...
                        %nufftWarp(pyr1(:, recR, c) .* AA1, ...
                        %        depChg12(:, recR), -depthClp1(:, recR), pind, filters, 2, 8) .* (1-mixFactor2(:, recR));
            end
        end

        pyrRec = pyrRec .* mulCorr;

        for iA = 1:appWrap,
            %sig = sigWrap * ones(nF, length(recR));
            %for k = nLev:-1:1,
            %    ind = pyrBandIndices(pind,k);
            %    sig(ind, :) = sig(ind, :) / (2^(nLev - k));%, appWrap/2);
            %end
            %pref = repmat(0.5 + 0.5 * erf((iA - 0.5) * dA./ (sqrt(2) * sig(:, 1))), [1, length(recR)]);

            pyrRec(:, :, appWrap + iA) = pyrRec(:, :, appWrap + iA)  + ...
                    pyrRec(:, :, 2*nA + appWrap + iA) ;
            pyrRec(:, :, 2*nA + appWrap + 1 - iA) = pyrRec(:, :, 2*nA + appWrap + 1 - iA)  + ...
                    pyrRec(:, :, appWrap + 1 - iA) ;
        end

        for iA = 1:2*nA,
            for iI = 1:length(recR),
                recons = reconSCFpyrGen1D(pyrRec(:, iI, iA+appWrap), pind, filters, 1, nLev) ;
                synth(iI, :, c, iA) = recons;
            end
        end
    end
    %synth(:, :, :, 1) = abs(synth(:, :, :, 1) - img1) * 4;
    %synth(:, :, :, 2) = abs(synth(:, :, :, 2) - img2) * 4;
    if false,
        %XXX resize back
        for iA = 1:2*nA,
            synth(:, :, :, iA) = imresize(synthT(:, :, :, iA), [nR, nCf]);
        end
    end
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
    %if lev == 1,
    %    curFilter = 2 * circshift(curFilter, [0, -ceil(lenO / 2)+1]);
    %    ctr = ceil((n+0.5)/2);
    %    curFilter(ctr) = curFilter(ctr)/2;
    %elseif lev ~= nLev,
    %    curFilter = 2 * circshift(curFilter, [0, -ceil(lenO / 4)]);
    %end

    curFilter = fftshift(fft(abs(ifft(ifftshift(curFilter)))));

    %plot(real(ifft(ifftshift(curFilter + [0, rot90(curFilter(2:n),2)])))); %XXX
    %plot(real(curFilter + [0, rot90(curFilter(2:n),2)])); %XXX
    tempDFT = curFilter.*imdft; % Transform domain

    % shift to avoid cross-boundary wavelet
    tempDFT = tempDFT .* exp(2*pi*1i*((0:n-1) - floor(n/2))/nUnex);
    tempDFT = tempDFT(indices);
    curResult = ifft(ifftshift(tempDFT)) * length(indices) / nUnex;
    curResult = circshift(curResult, [0, floor(length(indices)/(2*nUnex))]);
    pyr = real(curResult(:));
end

function [resPyr, movVis] = nufftWarp(pyr, mov, depth, pind, filters, m, Q, depPhsRatios)
    %http://www.cscamm.umd.edu/programs/fam04/qing_liu_fam04.pdf
    nLev = max(size(filters));
    nR = size(pyr, 2);
    n = length(filters{1});
    resPyr = zeros(size(pyr));
    movVis.mov = cell(nLev, 1);
    movVis.sig = cell(nLev, 1);
    movVis.sigFilt = cell(nLev, 1);
    movVis.layer = cell(nLev, 1);
    movVis.recon = cell(nLev, 1);
    for k = nLev:-1:1,
        curFilter = filters{k};
        indices = pyrBandIndices(pind, k);

        %if k == nLev,
        %    mov(indices, :) = 0;
        %end
        if k == nLev,
            resPyr(indices, :) = pyr(indices, :);
            continue;
        end
        nL = pind(k, 2);
        q = floor(min(Q, nL-1)/2)*2;

        if false%k >= nLev - 5,
            resPyr(indices, :) = pyr(indices, :) .* exp(-i * mov(indices, :) / depPhsRatios{k});
            movVis.mov{k} = zeros(size(pyr(indices, 1)));
            movVis.sigFilt{k} = resPyr(indices, 1);
            movVis.recon{k} = reconSCFpyrGen1D(resPyr(:, 1), pind, filters, k, nLev);
            movVis.layer{k} = reconSCFpyrGen1D(resPyr(:, 1), pind, filters, k, k);
            continue;
        end

        %pyr(indices, :) = pyr(indices, :) .* (1 - atan((pi / 2)*(mov(circshift(indices(:), 1), :) - mov(circshift(indices(:), -1), :)) / (n / nL)) / (pi / 2));

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

        A = zeros(q+1, nL, nR);
        mt = m*(mov(indices, :)*nL/n + repmat((0:nL-1)', [1, nR]));
        mtFloor = round(mt);
        stack = zeros(nL, 1);
        W = ones(nL, nR);

        movVis.mov{k} = mod(mt+0.5*m, nL*m);
        if false,
            for iR = 1:nR,
                top = 1;
                for iL = 1:nL,
                    while top > 1 && mtFloor(iL, iR) <= mtFloor(stack(top-1), iR) ...
                            && depth(indices(iL), iR) > depth(indices(stack(top-1)), iR),
                        if k ~= nLev,
                            W(stack(top-1), iR) = 0;
                        end
                        top = top - 1;
                    end
                    if top == 1 || mtFloor(iL, iR) > mtFloor(stack(top-1), iR),
                        stack(top) = iL;
                        top = top + 1;
                    else,
                        if k ~= nLev,
                            W(iL, iR) = 0;
                        end
                    end
                end
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

        X = reshape(F * reshape(A, [q+1, nL*nR]), [q+1, nL, nR]);
        sig = zeros(m*nL+q, nR);
        sig2 = zeros(m*nL+q, nR);
        for iQ = 1:q+1,
            %idx = mtFloor+iQ + repmat((0:nR-1)*(m*nL+q), [nL, 1]);
            subs = [reshape(mtFloor+iQ, [nL*nR, 1]), reshape(repmat(1:nR, [nL, 1]), [nL*nR, 1])];
            sig = sig + accumarray(subs, reshape(pyr(indices, :) .* W .* reshape(X(iQ, :, :), [nL, nR]),...
                    [nL*nR, 1]), [m*nL+q, nR]);
            %sig(idx) = sig(idx) + pyr(indices, :) .* W .* squeeze(X(iQ, :, :));
            if false,%k == nLev,
                sig2 = sig2 + accumarray(subs, reshape(abs(W .* squeeze(X(iQ, :, :))),...
                        [nL*nR, 1]), [m*nL+q, nR]);

                %sig2(idx) = sig2(idx) + W .* squeeze(X(iQ, :, :));
            end
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
        movVis.sig{k} = sig(:, 1);
        if false,%k == nLev,
            sig2(1+m*nL:q/2+m*nL, :) = sig2(1+m*nL:q/2+m*nL, :) + sig2(1:q/2, :);
            sig2(q/2+1:q, :) = sig2(q/2+1:q, :) + sig2(q/2+1+m*nL:q+m*nL, :);
            sig2 = sig2(q/2+1:q/2+m*nL, :);
        end

        for iR = 1:nR,
            fftSig = fftshift(fft(sig(:, iR)));
            fftSig = fftSig ./ cos(pi*((0:m*nL-1)-m*nL/2)/(m*nL))';
            resPyr(indices, iR) = ifft(ifftshift(fftSig(m*nL/2-floor(nL/2)+(1:nL))));
            if k == nLev,
                if false,
                    fftSig2 = fftshift(fft(sig2(:, iR)));
                    fftSig2 = fftSig2 ./ cos(pi*((0:m*nL-1)-m*nL/2)/(m*nL))';
                    dem = ifft(ifftshift(fftSig2(m*nL/2-floor(nL/2)+(1:nL))));
                    dem(dem == 0) = 1;
                    resPyr(indices, iR) = resPyr(indices, iR) ./ sqrt(abs(dem));
                else
                    resPyr(indices, iR) = resPyr(indices, iR);
                end
            end
        end
        movVis.sigFilt{k} = resPyr(indices, 1);
        movVis.recon{k} = reconSCFpyrGen1D(resPyr(:, 1), pind, filters, k, nLev);
        movVis.layer{k} = reconSCFpyrGen1D(resPyr(:, 1), pind, filters, k, k);
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
