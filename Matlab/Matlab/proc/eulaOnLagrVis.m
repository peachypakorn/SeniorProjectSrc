function [synth, vis] = eulaOnLagrVis(img1, img2, depthLo, depthRo, testParam, runParam, noEula)
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
    if isfield(runParam, 'forcematch'),
        sigWrap = runParam.sigWrap;
    else,
        sigWrap = 0;
    end
    if isfield(runParam, 'recR'),
        recR = runParam.recR;
    else,
        recR = 0;
    end

    nR = size(img1, 1);
    nC = size(img1, 2);
    nCol = size(img1, 3);
    synth = zeros([ nR, nC, nCol,  nA*2]);
    maxAmp = (nA-1)*dA;
    sigma_v_fac = 1;

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

    pindRef = buildBandIndices(filters, 4);
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


        levDepthL = imresize(lagrDepthL, [nR, pind(iLev, 2)*4])';
        sigma_v = 2^(iLev-2) * sigma_v_fac;
        h = fspecial('gaussian', [1, ceil(sigma_v*2*2+1)], sigma_v+1e-10);
        levDepthL = imfilter(levDepthL, h, 'symmetric');
        levDepthL = reduce(@min, levDepthL, 4, 2);
        preDepthCont1(ind, :) = levDepthL/dsrLagr;

        preDepth1(ind, :) = repmat((-0.5 + (1:nO)')*(nC/nO), [1, nR]) - levDepthL/dsrLagr;
        preDepth1(ind, :) = clip(round(preDepth1(ind, :) * (nRef / nC) - refStart), 0, nRef-1);
        refInd1(ind, :, :) = repmat(preDepth1(ind, :) + pindRef(iLev, 1), [1, 1, nCol])+ ...
                             repmat((0:nR-1)*nFRef, [length(ind), 1, nCol]) + ...
                             repmat(reshape((0:nCol-1)*nFRef*nR, [1, 1, nCol]), [length(ind), nR, 1]);
        preDepth1(ind, :) = (preDepth1(ind, :) + refStart) * (nC / nRef) - ...
                repmat((-0.5 + (1:nO)')*(nC/nO), [1, nR]);

        levDepthR = imresize(lagrDepthR, [nR, pind(iLev, 2)*4])';
        levDepthR = imfilter(levDepthR, h, 'symmetric');
        levDepthR = reduce(@min, levDepthR, 4, 2);
        preDepthCont2(ind, :) = levDepthR/dsrLagr;

        preDepth2(ind, :) = repmat((-0.5 + (1:nO)')*(nC/nO), [1, nR]) + levDepthR/dsrLagr;
        preDepth2(ind, :) = clip(round(preDepth2(ind, :) * (nRef / nC) - refStart), 0, nRef-1);
        refInd2(ind, :, :) = repmat(preDepth2(ind, :) + pindRef(iLev, 1), [1, 1, nCol]) + ...
                             repmat((0:nR-1)*nFRef, [length(ind), 1, nCol]) + ...
                             repmat(reshape((0:nCol-1)*nFRef*nR, [1, 1, nCol]), [length(ind), nR, 1]);
        preDepth2(ind, :, :) = repmat((-0.5 + (1:nO)')*(nC/nO), [1, nR]) - ...
                (preDepth2(ind, :) + refStart) * (nC / nRef);
    end

    for iR=1:nR,
        for iCol = 1:nCol,
            pyr1(:, iR, iCol) = buildSCFpyrGen1D(img1(iR, :, iCol), filters, 1);
            pyr1Ref(:, iR, iCol) = buildSCFpyrGen1D(img1(iR, :, iCol), filters, 4);
            pyr2(:, iR, iCol) = buildSCFpyrGen1D(img2(iR, :, iCol), filters, 1);
            pyr2Ref(:, iR, iCol) = buildSCFpyrGen1D(img2(iR, :, iCol), filters, 4);
        end
    end;

    vis.pyr1 = cell(nLev, 1);
    vis.pyr2 = cell(nLev, 1);
    vis.layer1 = cell(nLev, 1);
    vis.layer2 = cell(nLev, 1);
    vis.accum1 = cell(nLev, 1);
    vis.accum2 = cell(nLev, 1);
    accum1 = zeros(1, nC);
    accum2 = zeros(1, nC);
    for iLev = nLev:-1:1,
        ind = pyrBandIndices(pind, iLev);
        vis.pyr1{iLev} = pyr1(ind, 1, 1);
        vis.pyr2{iLev} = pyr2(ind, 1, 1);

        tmp = zeros(nF, 1);
        tmp(ind) = pyr1(ind, 1, 1);
        tmp = reconSCFpyrGen1D(tmp, pind, filters);
        vis.layer1{iLev} = tmp;
        accum1 = accum1 + tmp;
        vis.accum1{iLev} = accum1;

        tmp = zeros(nF, 1);
        tmp(ind) = pyr2(ind, 1, 1);
        tmp = reconSCFpyrGen1D(tmp, pind, filters);
        vis.layer2{iLev} = tmp;
        accum2 = accum2 + tmp;
        vis.accum2{iLev} = accum2;
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

    sigma_v = 2^(nLev-2) * sigma_v_fac;
    ind = pyrBandIndices(pind, nLev);
    h = fspecial('gaussian', [ceil(sigma_v*2*2+1), 1], sigma_v+1e-10);
    tmpDepthL = -imfilter(imresize(lagrDepthL, [nR, pind(nLev, 2)*4]), h, 'symmetric')'/dsrLagr;
    preDepth1(ind, :) = reduce(@max, tmpDepthL, 4, 2);
    tmpDepthR = -imfilter(imresize(lagrDepthR, [nR, pind(nLev, 2)*4]), h, 'symmetric')'/dsrLagr;
    preDepth2(ind, :) = reduce(@max, tmpDepthR, 4, 2);
    phsDif1(ind, :, :) = 0;
    phsDif2(ind, :, :) = 0;
    weight1(ind, :, :) = 1;
    weight2(ind, :, :) = 1;

    for k = nLev:-1:1,
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
        sigma_v = cLen * sigma_v_fac + 1e-10;
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

    %figure;
    %hold on;
    %ind = pyrBandIndices(pind,8);
    %plot(phsDif1(ind, 10, 1));
    %plot(phsDif1(ind, 9, 1));
    %plot(phsDif1(ind, 11, 1));
    %plot(preDepth1(ind, 10)/ depPhsRatios{8}, 'black');
    %plot((1:nC*dsrLagr) * pind(8, 2) / (nC*dsrLagr), -lagrDepthL(10, :)'/(dsrLagr*depPhsRatios{8}), 'red');

    %figure;
    %hold on;
    %plot(abs(pyr1(ind, 10, 1)));
    %plot(abs(pyr2(ind, 10, 1)), 'red');
    %plot(abs(pyr2Ref(refInd1(ind, 10, 1))), 'black');


    %error('1324')
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
    if false,
    %between-layer prop
        for k = nLev-2:-1:1,
            ind = pyrBandIndices(pind,k);
            ind_l = pyrBandIndices(pind,k+1);
            for iR = 1:nR,
                for idx = 1:max(size(ind_l)),
                    if k > 1,
                        for d = -1:0,
                            if idx*2+d > 0 && idx*2+d <= max(size(ind)),
                                if abs(phsDif1(ind_l(idx), iR)) >= pi / 4,
                                    phsDif1(ind(idx*2+d), iR) = phsDif1(ind_l(idx), iR) * 2;
                                    mix = clip((abs(phsDif1(ind_l(idx), iR)) - pi / 4) / (pi / 4), 0, 1);
                                    depth1(ind(idx*2+d), iR) = depth1(ind_l(idx), iR) * mix + depth1(ind(idx*2+d), iR) * (1-mix);
                                end
                                if abs(phsDif2(ind_l(idx), iR)) >= pi / 4,
                                    phsDif2(ind(idx*2+d), iR) = phsDif2(ind_l(idx), iR) * 2;
                                    mix = clip((abs(phsDif2(ind_l(idx), iR)) - pi / 4) / (pi / 4), 0, 1);
                                    depth2(ind(idx*2+d), iR) = depth2(ind_l(idx), iR) * mix + depth2(ind(idx*2+d), iR) * (1-mix);
                                end
                            end
                        end
                    else,
                        if abs(phsDif1(ind_l(idx), iR)) >= pi / 4,
                            mix = clip((abs(phsDif1(ind_l(idx), iR)) - pi / 4) / (pi / 4), 0, 1);
                            depth1(ind(idx), iR) = depth1(ind_l(idx), iR) * mix + depth1(ind(idx), iR) * (1-mix);
                        end
                        if abs(phsDif2(ind_l(idx), iR)) >= pi / 4,
                            mix = clip((abs(phsDif2(ind_l(idx), iR)) - pi / 4) / (pi / 4), 0, 1);
                            depth2(ind(idx), iR) = depth2(ind_l(idx), iR) * mix + depth2(ind(idx), iR) * (1-mix);
                        end
                    end
                end
            end
        end
    end

    if maxDepth > 0,
        depthClp1 = maxDepth * atan(depth1 / maxDepth);
        depthClp2 = maxDepth * atan(depth2 / maxDepth);
        %depthClp = clip(depth, -maxDepth, maxDepth);
    else
        depthClp1 = depth1;
        depthClp2 = depth2;
    end

    if false,
    %forcematch
        coseDif1 = zeros(nF, nR);
        coseDif2 = zeros(nF, nR);
        weight1 = zeros(nF, nR);
        weight2 = zeros(nF, nR);
        for c = 1:nCol,
            for iLev = nLev-1:-1:1,
                ind = pyrBandIndices(pind, iLev);
                tv = exp(1i * (angle(pyr1(2:pind(iLev, 2), :, c)) - angle(pyr1(1:pind(iLev, 2)-1, :, c))));
                tw = bsxfun(@min, abs(pyr1(2:pind(iLev, 2), :, c)), abs(pyr1(1:pind(iLev, 2)-1, :, c)));
                coseDif1(ind(2:pind(iLev, 2)), :) = coseDif1(ind(2:pind(iLev, 2)), :) + tv .* tw;
                weight1(ind(2:pind(iLev, 2)), :) = weight1(ind(2:pind(iLev, 2)), :) + tw;
                coseDif1(ind(1:pind(iLev, 2)-1), :) = coseDif1(ind(1:pind(iLev, 2)-1), :) + tv .* tw;
                weight1(ind(1:pind(iLev, 2)-1), :) = weight1(ind(1:pind(iLev, 2)-1), :) + tw;

                tv = exp(1i * (angle(pyr2(2:pind(iLev, 2), :, c)) - angle(pyr2(1:pind(iLev, 2)-1, :, c))));
                tw = bsxfun(@min, abs(pyr2(2:pind(iLev, 2), :, c)), abs(pyr2(1:pind(iLev, 2)-1, :, c)));
                coseDif2(ind(2:pind(iLev, 2)), :) = coseDif2(ind(2:pind(iLev, 2)), :) + tv .* tw;
                weight2(ind(2:pind(iLev, 2)), :) = weight2(ind(2:pind(iLev, 2)), :) + tw;
                coseDif2(ind(1:pind(iLev, 2)-1), :) = coseDif2(ind(1:pind(iLev, 2)-1), :) + tv .* tw;
                weight2(ind(1:pind(iLev, 2)-1), :) = weight2(ind(1:pind(iLev, 2)-1), :) + tw;
            end
        end
        coseDif1 = angle(coseDif1 ./ (weight1 + 1e-10));
        coseDif2 = angle(coseDif2 ./ (weight2 + 1e-10));
        matDepth1 = 1e10 * ones([nF, nR]);
        matDepth2 = 1e10 * ones([nF, nR]);
        for iLev = nLev-1:-1:1,
            ind = pyrBandIndices(pind, iLev);
            matDepth1(ind, :) = (2*pi ./ (abs(coseDif1(ind, :)) + 1e-10)) * nC / length(ind);
            matDepth2(ind, :) = (2*pi ./ (abs(coseDif2(ind, :)) + 1e-10)) * nC / length(ind);
        end
        fac = 2*nA*dA;
        valDepth1 = round(depthClp1 * fac ./ matDepth1) .* matDepth1 / fac;
        valDepth2 = round(depthClp2 * fac ./ matDepth2) .* matDepth2 / fac;
        %depthClp1 = valDepth1;
        %depthClp2 = valDepth2;
        for iLev = nLev-1:-1:1,
            ind = pyrBandIndices(pind, iLev);
            %ind
            %sum(abs(depthClp1(ind, :))  * fac > matDepth1(ind, :))
        end
        depthClp1(abs(depthClp1)  * fac > matDepth1/2) = valDepth1(abs(depthClp1)  * fac > matDepth1/2);
        depthClp2(abs(depthClp2) * fac > matDepth2/2) = valDepth2(abs(depthClp2)  * fac > matDepth2/2);
    end
    %for k = nLev-1:-1:1,
    %    ind = pyrBandIndices(pind,k);
    %    indChg = abs(depthClp(ind) / depPhsRatios{k}) > pi/(2*nA*dA);
    %    depthClp(ind(indChg)) = round(depthClp(ind(indChg)) * 2*nA*dA / ...
    %            (depPhsRatios{k} * 2 * pi)) * (pi/(nA*dA)) * depPhsRatios{k};
    %end

    phsDifAA1 = zeros([nF, nR]);
    phsDifAA2 = zeros([nF, nR]);
    for k = nLev-1:-1:1,
        ind = pyrBandIndices(pind,k);
        phsDifAA1(ind, :) = dA * depthClp1(ind, :) / depPhsRatios{k};
        phsDifAA2(ind, :) = dA * depthClp2(ind, :) / depPhsRatios{k};
    end

    if recR == 0,
        recR = 1:nR;
    end
    vis.mov = cell(2*nA, 2);
    for c = 1:nCol,
        %% Reconstructing
        appWrap = ceil(sigWrap*2);
        pyrRec = zeros(nF, nR, 2*nA+2*appWrap);
        for iA = 1:nA+appWrap,
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
            if c == 1,
                [tmp, movVis] = nufftWarp(pyr1(:, 1:1, c) .* exp( -( phsDifAA1(:, 1:1) * rho).^2 / 2), ...
                            -depChg11(:, 1:1), -depthClp1(:, 1:1), pind, filters, 2, 8);
                vis.mov{nA+1-iA, 1} = movVis;
                [tmp, movVis] = nufftWarp(pyr2(:, 1:1, c) .* exp( -( phsDifAA2(:, 1:1) * rho).^2 / 2), ...
                            -depChg22(:, 1:1), -depthClp2(:, 1:1), pind, filters, 2, 8);
                vis.mov{nA+1-iA, 2} = movVis;
                [tmp, movVis] = nufftWarp(pyr2(:, 1:1, c) .* exp( -( phsDifAA2(:, 1:1) * rho).^2 / 2), ...
                            depChg21(:, 1:1), -depthClp2(:, 1:1), pind, filters, 2, 8);
                vis.mov{nA+iA, 1} = movVis;
                [tmp, movVis] = nufftWarp(pyr1(:, 1:1, c) .* exp( -( phsDifAA1(:, 1:1) * rho).^2 / 2), ...
                            depChg12(:, 1:1), -depthClp1(:, 1:1), pind, filters, 2, 8);
                vis.mov{nA+iA, 2} = movVis;
            end

            pyrRec(:, recR, nA+appWrap+1-iA) = ...
                    nufftWarp(pyr1(:, recR, c) .* exp( -( phsDifAA1(:, recR) * rho).^2 / 2), ...
                            -depChg11(:, recR), -depthClp1(:, recR), pind, filters, 2, 8) .* mixFactor1(:, recR) +...
                    nufftWarp(pyr2(:, recR, c) .* exp( -( phsDifAA2(:, recR) * rho).^2 / 2), ...
                            -depChg22(:, recR), -depthClp2(:, recR), pind, filters, 2, 8) .* (1-mixFactor1(:, recR));

            pyrRec(:, recR, nA+appWrap+iA) = ...
                    nufftWarp(pyr2(:, recR, c) .* exp( -( phsDifAA2(:, recR) * rho).^2 / 2), ...
                            depChg21(:, recR), -depthClp2(:, recR), pind, filters, 2, 8) .* mixFactor2(:, recR) +...
                    nufftWarp(pyr1(:, recR, c) .* exp( -( phsDifAA1(:, recR) * rho).^2 / 2), ...
                            depChg12(:, recR), -depthClp1(:, recR), pind, filters, 2, 8) .* (1-mixFactor2(:, recR));
        end

        for iA = 1:appWrap,
            mixFactor = exp((iA - 0.5) / sigWrap) / (1 + exp((iA - 0.5) / sigWrap));
            pyrRec(:, :, appWrap + iA) = pyrRec(:, :, appWrap + iA) * mixFactor + ...
                    pyrRec(:, :, 2*nA + appWrap + iA) * (1-mixFactor);
            pyrRec(:, :, 2*nA + appWrap + 1 - iA) = pyrRec(:, :, 2*nA + appWrap + 1 - iA) * mixFactor + ...
                    pyrRec(:, :, appWrap + 1 - iA) * (1-mixFactor);
        end

        for iA = 1:2*nA,
            for iR = recR,
                recons = reconSCFpyrGen1D(pyrRec(:, iR, iA+appWrap), pind, filters) ;
                synth(iR, :, c, iA) = recons;
            end
        end
    end
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
        reduceArr(2:nR/fac, :) = bsxfun(fcn, reduceArr(2:nR/fac, :), arr(ind+4:fac:nR-4, :));
    end
    for ind = fac+1:fac+ext,
        reduceArr(1:nR/fac-1, :) = bsxfun(fcn, reduceArr(1:nR/fac-1, :), arr(ind:fac:nR, :));
    end
end

function outFilters = get1DRadialFilters(dims, rVals, twidth)
    % Construct log_rad
    ctr = ceil((dims+0.5)/2);
    log_rad = abs(((1:dims)-ctr)./(dims/2));
    log_rad(ctr) = log_rad(ctr-1);

    N = max(size(rVals));

    [himask, lomaskPrev] = getRadialMaskPair(rVals(1),log_rad, twidth);
    mask = (1:dims)>ctr;
    outFilters{1} = himask.*mask;
    for k = 2:N
        [himask, lomask] = getRadialMaskPair(rVals(k),log_rad, twidth);
        outFilters{k} = himask.*lomaskPrev.*mask;
        lomaskPrev = lomaskPrev.*lomask;
    end
    outFilters{N+1} = lomaskPrev;
end

function [himask, lomask] = getRadialMaskPair( r, log_rad, twidth)
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

function [resPyr, movVis] = nufftWarp(pyr, mov, depth, pind, filters, m, Q)
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
        nL = pind(k, 2);
        q = floor(min(Q, nL-1)/2)*2;
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
        if k ~= nLev,
            movVis.mov{k} = mod(mt+0.5*m, nL*m);
        else
            mt = m * repmat((0:nL-1)', [1, nR]);
            movVis.mov{k} = (0.5 + (0:nL-1)) * m;
        end
        mtFloor = round(mt);
        stack = zeros(nL, 1);
        W = ones(nL, nR);
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
                sig2 = sig2 + accumarray(subs, reshape(abs(W .* reshape(X(iQ, :, :), [nL, nR])),...
                        [nL*nR, 1]), [m*nL+q, nR]);

                %sig2(idx) = sig2(idx) + W .* squeeze(X(iQ, :, :));
            end
        end
        sig(1+m*nL:q/2+m*nL, :) = sig(1+m*nL:q/2+m*nL, :) + sig(1:q/2, :);
        sig(q/2+1:q, :) = sig(q/2+1:q, :) + sig(q/2+1+m*nL:q+m*nL, :);
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
                    resPyr(indices, iR) = pyr(indices, iR);
                end
            end
        end
        movVis.sigFilt{k} = resPyr(indices, 1);
        movVis.recon{k} = reconSCFpyrGen1D(resPyr(:, 1), pind, filters);
        tmp =zeros(size(pyr, 1), 1);
        tmp(indices, 1) = resPyr(indices, 1);
        movVis.layer{k} = reconSCFpyrGen1D(tmp, pind, filters);
    end
end

function res = reconSCFpyrGen1D(pyr, pind, filters)
    nLev = max(size(filters));
    n = length(filters{1});
    imdft = zeros(1, n); %DFT of image
    for k = 1:nLev
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

function [ out ] = getIDXFromFilter1D(filt, expand)
    aboveZero = filt>1e-10;
    aboveZero = logical(aboveZero + fliplr(aboveZero));
    n = length(filt);

    idx1 = 1:length(filt);
    idx1 = idx1(aboveZero);
    len = max(idx1) - min(idx1) + 1;

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
