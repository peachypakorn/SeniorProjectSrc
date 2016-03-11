function [synth, mem, phsCmpO] = eulaOnLagrGuideProp(img1, img2, depthLo, depthRo, testParam, runParam, mem, iA, noEula, phsCmpO)
    rho = runParam.rho;
    ht = runParam.ht;
    nA = runParam.nA;
    dA = runParam.dA;
    dsrLagr = runParam.eurLagrResize;
    maxDLagr = ceil(testParam.maxDepth * dsrLagr);
    vis = struct();
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

    lagrDepthL = depthLo;
    lagrDepthR = depthRo;
    dsrLagr = size(lagrDepthL, 2) / size(img1, 2);

    waveletExt = 4;

    nR = size(img1, 1);
    nC = size(img1, 2);
    nCol = size(img1, 3);
    maxAmp = (nA-1)*dA;
    sigma_v_fac = 1e-6;

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
                pyr2(:, iR, iCol) = buildSCFpyrGen1D(img2(iR, :, iCol), filters, 1);

                pyr1Ref(:, iR, iCol) = buildSCFpyrGen1D(img1(iR, :, iCol), filters, waveletExt);
                pyr2Ref(:, iR, iCol) = buildSCFpyrGen1D(img2(iR, :, iCol), filters, waveletExt);
            end
        end;

        depth1 = zeros([nF, nR]);
        depth2 = zeros([nF, nR]);
        levDepthL = [0];
        levDepthR = [0];
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
            if size(lagrDepthL, 2)>1,
                levDepthLU = imresize(levDepthL, [nO*4, nR]);
                levDepthRU = imresize(levDepthR, [nO*4, nR]);

                levDepthL = imresize(lagrDepthL, [nR, pind(iLev, 2)*4])'/dsrLagr;
                mask = abs(levDepthL - levDepthLU) < (1.0 /  dsrLagr);
                fact = max(0, abs(levDepthL - levDepthLU) * dsrLagr * 2 - 1);
                levDepthL(mask) = levDepthLU(mask) .* (1 - fact(mask)) + levDepthL(mask) .* fact(mask);
                levDepthL = reduce(@min, levDepthL, 4, 0);

                levDepthR = imresize(lagrDepthR, [nR, pind(iLev, 2)*4])'/dsrLagr;
                mask = abs(levDepthR - levDepthRU) < (1.0 /  dsrLagr);
                fact = max(0, abs(levDepthR - levDepthRU) * dsrLagr * 2 - 1);
                levDepthR(mask) = levDepthRU(mask) .* (1 - fact(mask)) + levDepthR(mask) .* fact(mask);
                levDepthR = reduce(@min, levDepthR, 4, 0);
            else,
                levDepthL = imresize(levDepthL, [nO, nR]);
                levDepthR = imresize(levDepthR, [nO, nR]);
            end

            preDepth1 = repmat((-0.5 + (1:nO)')*(nC/nO), [1, nR]) - levDepthL;
            preDepth1 = clip(round(preDepth1 * (nRef / nC) - refStart), 0, nRef-1);
            refInd1 = repmat(preDepth1 + pindRef(iLev, 1), [1, 1, nCol]);
            preDepth1 = (preDepth1 + refStart) * (nC / nRef) - ...
                    repmat((-0.5 + (1:nO)')*(nC/nO), [1, nR]);

            preDepth2 = repmat((-0.5 + (1:nO)')*(nC/nO), [1, nR]) + levDepthR;
            preDepth2 = clip(round(preDepth2 * (nRef / nC) - refStart), 0, nRef-1);
            refInd2 = repmat(preDepth2 + pindRef(iLev, 1), [1, 1, nCol]);
            preDepth2 = repmat((-0.5 + (1:nO)')*(nC/nO), [1, nR]) - ...
                    (preDepth2 + refStart) * (nC / nRef);

            %% Computing Filter
            phsDif1 = zeros([nO nR nCol]);
            phsDif2 = zeros([nO nR nCol]);
            weight1 = zeros([nO nR nCol]);
            weight2 = zeros([nO nR nCol]);
            for iR=1:nR,
                for iCol = 1:nCol,
                    phsDif1(:, iR, iCol) = angle(pyr2Ref(refInd1(:, iR, iCol), iR, iCol)) - angle(pyr1(ind, iR, iCol));
                    phsDif2(:, iR, iCol) = angle(pyr2(ind, iR, iCol)) - angle(pyr1Ref(refInd2(:, iR, iCol), iR, iCol));
                    weight1(:, iR, iCol) = bsxfun(@min, abs(pyr2Ref(refInd1(:, iR, iCol), iR, iCol)), abs(pyr1(ind, iR, iCol)));
                    weight2(:, iR, iCol) = bsxfun(@min, abs(pyr2(ind, iR, iCol)), abs(pyr1Ref(refInd2(:, iR, iCol), iR, iCol)));
                end
            end
            phsDif1 = angle(sum(exp(phsDif1*1i) .* weight1(:, :, :), 3));
            phsDif2 = angle(sum(exp(phsDif2*1i) .* weight2(:, :, :), 3));
            weight1 = sum(weight1, 3);
            weight2 = sum(weight2, 3);

            levDepthL = -preDepth1 + phsDif1 * depPhsRatios{iLev};
            levDepthR = -preDepth2 + phsDif2 * depPhsRatios{iLev};

            sigma_h = max(0, sqrt(12 / (nR / nO)));
            sigma_v = sigma_h * nR / nO;
            kH = ceil(sigma_h*4);
            kV = ceil(sigma_v*4);

            h_h = fspecial('gaussian', [kH, 1], sigma_h);
            h_v = fspecial('gaussian', [1, kV], sigma_v);
            sumDep1 = imfilter(imfilter(levDepthL .* weight1, h_h, 'circular'), h_v, 'symmetric');
            sumDep2 = imfilter(imfilter(levDepthR .* weight2, h_h, 'circular'), h_v, 'symmetric');
            wDep1 = imfilter(imfilter(weight1, h_h, 'circular'), h_v, 'symmetric')+1e-20;
            wDep2 = imfilter(imfilter(weight2, h_h, 'circular'), h_v, 'symmetric')+1e-20;

            levDepthL = sumDep1 ./ wDep1;
            levDepthR = sumDep2 ./ wDep2;
            depth1(ind, :) = levDepthL;
            depth2(ind, :) = levDepthR;
            imshow(imresize(imresize(levDepthL' / 100, [nR, nC]), 1) + 0.5)
            figure
        end

        if maxDepth > 0,
            depthClp1 = maxDepth * atan(depth1 / maxDepth);
            depthClp2 = maxDepth * atan(depth2 / maxDepth);
        else
            depthClp1 = depth1;
            depthClp2 = depth2;
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

        amp = (nA+0.5-iA)*dA;

        depChg1 = mem.depthClp1 * amp - mem.depth1/2;
        depChg2 = mem.depthClp2 * amp + mem.depth2/2;
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


        if expSinc,
            AA1 = mem.phsDifAA1(:, recR) * rho < pi/2;
            AA2 = mem.phsDifAA2(:, recR) * rho < pi/2;
        else,
            AA1 = exp( -( mem.phsDifAA1(:, recR) * rho).^2 / 2);
            AA2 = exp( -( mem.phsDifAA2(:, recR) * rho).^2 / 2);
        end

        %if iA <= nA,
        if amp > 0.5,
            pyrRec(:, :) = ...
                    nufftWarp(mem.pyr1(:, recR, c) .* AA1, ...
                    depChg1(:, recR), -mem.depthClp1(:, recR), pind, filters, 2, 4, depPhsRatios, -1);
        elseif amp > -0.5,
            [synth1, weight1] = ...
                    nufftWarp(mem.pyr1(:, recR, c) .* AA1, ...
                    depChg1(:, recR), -mem.depthClp1(:, recR), pind, filters, 2, 4, depPhsRatios, 1);
            [synth2, weight2] = ...
                    nufftWarp(mem.pyr2(:, recR, c) .* AA2, ...
                    depChg2(:, recR), -mem.depthClp2(:, recR), pind, filters, 2, 4, depPhsRatios, -1);
            pyrRec(:, :) = (synth1 .* weight1 + synth2 .* weight2) ./ (weight1 + weight2);
        else
            pyrRec(:, :) = ...
                    nufftWarp(mem.pyr2(:, recR, c) .* AA2, ...
                    depChg2(:, recR), -mem.depthClp2(:, recR), pind, filters, 2, 4, depPhsRatios, 1);
        end
        
        pyrRec = pyrRec;

        for iI = 1:length(recR),
            recons = reconSCFpyrGen1D(pyrRec(:, iI), pind, filters, 1, nLev) ;
            synth(iI, :, c) = recons;
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

    curFilter = fftshift(fft(abs(ifft(ifftshift(curFilter)))));

    tempDFT = curFilter.*imdft; % Transform domain

    % shift to avoid cross-boundary wavelet
    tempDFT = tempDFT .* exp(2*pi*1i*((0:n-1) - floor(n/2))/nUnex);
    tempDFT = tempDFT(indices);
    curResult = ifft(ifftshift(tempDFT)) * length(indices) / nUnex;
    curResult = circshift(curResult, [0, floor(length(indices)/(2*nUnex))]);
    pyr = real(curResult(:));
end

function [resPyr, weight] = nufftWarp(pyr, mov, depth, pind, filters, m, Q, depPhsRatios, direction)
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
        for k = nLev:-1:1,
            curFilter = filters{k};
            indices = pyrBandIndices(pind, k);

            %if k == nLev,
            %    mov(indices, :) = 0;
            %end
            if k == nLev,
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

            if direction < 0,
                cum = zeros(nL, rSize);
                curMin = 1e9*ones(1, rSize);
                for iL = nL:-1:1,
                    cum(iL, :) = curMin;
                    curMin = min(curMin, mt(iL, :));
                end
                %fac = exp(-2*(max((mt - cum)/m + 1, 0)));
                fac = cos(pi * (min(max((mt - cum)/m + 1, 0), 0.5)));
                pyr(indices, rInd) = pyr(indices, rInd) .* fac;
            else,
                cum = zeros(nL, rSize);
                curMax = -1e9*ones(1, rSize);
                for iL = 1:nL,
                    cum(iL, :) = curMax;
                    curMax = max(curMax, mt(iL, :));
                end
                %fac = exp(-2*(max((cum - mt)/m + 1, 0)));
                fac = cos(pi * (min(max((cum - mt)/m + 1, 0), 0.5)));
                pyr(indices, rInd) = pyr(indices, rInd) .* fac;
            end

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
                w2 = (cos((pi / 4) * min(2, max(0, dif - 1))) + 1e-3);
                wr = interp1(movs, [w1, w2], (0:nL-1)', 'pchip');
                weight(indices, rId) = wr(:, 2) ./ wr(:, 1);
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

function filtered = bilateral(img, guid, sig_h, sig_v, sig_c)
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
    frameEx(vEx + (1:nR), hEx + (1:nC)) = 1;

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
