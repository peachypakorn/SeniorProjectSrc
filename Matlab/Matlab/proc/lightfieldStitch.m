function [phsCmpO] = lightfieldStitch(ifiles, oprefix, runParam)
    ht = runParam.ht;
    nA = length(ifiles);

    sigWrap = runParam.sigWrap;
    appWrap = runParam.appWrap;
    reduWrap = runParam.reduWrap;
    phsCmpO = 0;

    if isfield(runParam, 'expShear'),
        expShear = runParam.expShear;
    else,
        expShear = false;
    end
    keyLev = runParam.shearLev;

    phsCmp = 0;

    startA = 1+appWrap;
    corrDist = nA - appWrap + 1 - startA;
    img1 = double(imread(ifiles{startA})) / 255;
    img2 = double(imread(ifiles{min(startA+corrDist, nA)})) / 255;

    nR = size(img1, 1);
    nC = size(img1, 2);
    nCol = size(img1, 3);

    sizFilt = 2.^[0:-1:-ht];
    filters = get1DRadialFilters(nC, sizFilt, 1);
    nLev = max(size(filters));
    depPhsRatios = getDepthPhaseRatio(filters);

    pind = buildBandIndices(filters, 1);
    nF = pind(nLev, 1) + pind(nLev, 2) - 1;

    pyrCmp1 = zeros([nF nR nCol]);
    pyrCmp2 = zeros([nF nR nCol]);

    for iCol = 1:nCol,
        for iR=1:nR,
            pyrCmp1(:, iR, iCol) = buildSCFpyrGen1D(img1(iR, :, iCol), filters, 1, runParam.shearMinLev, runParam.shearLev);
            pyrCmp2(:, iR, iCol) = buildSCFpyrGen1D(img2(iR, :, iCol), filters, 1, runParam.shearMinLev, runParam.shearLev);
        end;
    end
    for iLev = runParam.shearMinLev:runParam.shearLev,
        ind = pyrBandIndices(pind, iLev);
        if false
            fac = bsxfun(@max, abs(pyrCmp1(ind, :, :)), abs(pyrCmp2(ind, :, :))) ./ ...
                  bsxfun(@min, abs(pyrCmp1(ind, :, :)), abs(pyrCmp2(ind, :, :)));
            pyrCmp1(ind, :, :) = pyrCmp1(ind, :, :) ./ fac;
            pyrCmp2(ind, :, :) = pyrCmp2(ind, :, :) ./ fac;
        end
    end
    pyrCmp1 = pyrCmp1(:, :, 1) * 0.3 + pyrCmp1(:, :, 2) * 0.6 + pyrCmp1(:, :, 3) * 0.1;
    pyrCmp2 = pyrCmp2(:, :, 1) * 0.3 + pyrCmp2(:, :, 2) * 0.6 + pyrCmp2(:, :, 3) * 0.1;

    if expShear,
        phsCmpA = 0;
        for curLev = keyLev:-1:runParam.shearMinLev,
            po2 = 2^(keyLev - curLev);
            gR = min(runParam.shearGrid(1) * po2, nR);
            gC = runParam.shearGrid(2) * po2;
            inds = pyrBandIndices(pind, curLev);
            nL = length(inds);
            phsCont = (nC / nL) / depPhsRatios{curLev};
            imCmp1 = pyrCmp1(inds, :) .* repmat(exp(-i*(1:nL)'*phsCont), [1, nR]);
            imCmp2 = pyrCmp2(inds, :) .* repmat(exp(-i*(1:nL)'*phsCont), [1, nR]);
            imCmp1 = imresize(imCmp1, [gC, nR]);
            imCmp2 = imresize(imCmp2, [gC, nR]);
            dif = imCmp1 ./ imCmp2;
            dif = dif .* bsxfun(@min, abs(imCmp1), abs(imCmp2)) ./ abs(dif);
            phsCmp = angle(imresize(dif, [gC, gR])) / (corrDist * po2);
            if ~isscalar(phsCmpA),
                phsCmpA = imresize(phsCmpA, [gC, gR]);
                phsCmpA = phsCmp + (2*pi / (corrDist*po2)) * round((phsCmpA - phsCmp) / (2*pi / (corrDist*po2)));
            else
                phsCmpA = phsCmp;
            end
        end
        phsCmpA = -phsCmpA;
        phsCmp = phsCmpA;
        phsCmpO = phsCmpA;
        mkdir(sprintf('%s/sup/', oprefix));
        imwrite(imresize(0.5 + phsCmpO'*corrDist/(4*pi), [nR, nC]), sprintf('%s/sup/phsCmp.jpg', oprefix),'jpg');
    end

    if appWrap-reduWrap > 0,
        inds = pyrBandIndices(pind, keyLev+1);
        nCE = length(inds);
        nRE = ceil(nR * nCE / nC);
        errMap = zeros(2*(appWrap - reduWrap), nRE, nCE);
        for iA = 1+reduWrap:2*appWrap-reduWrap,
            iA
            img1 = double(imread(ifiles{iA})) / 255;
            pyr1 = zeros([nF, nR, nCol]);

            for iR=1:nR,
                for iCol = 1:nCol,
                    pyr1(:, iR, iCol) = buildSCFpyrGen1D(img1(iR, :, iCol), filters, 1, keyLev+1, keyLev+1);
                end
            end

            if expShear,
                disCorr = zeros(nF, nR);
                for iLev = 1:nLev-1,
                    inds = pyrBandIndices(pind, iLev);
                    imCorr = phsCmp * (iA-0.5*(startA+nA-appWrap)) * depPhsRatios{keyLev};
                    disCorr(inds, :) = imresize(imCorr, [length(inds), nR]);
                end
                for iCol = 1:nCol,
                    pyr1(:, :, iCol) = nufftWarp(pyr1(:, :, iCol), -disCorr, pind, filters, 2, 8, keyLev+1, keyLev+1);
                end
            end

            iA2 = iA + corrDist;
            img2 = double(imread(ifiles{iA2})) / 255;
            pyr2 = zeros([nF, nR, nCol]);

            for iR=1:nR,
                for iCol = 1:nCol,
                    pyr2(:, iR, iCol) = buildSCFpyrGen1D(img2(iR, :, iCol), filters, 1, keyLev+1, keyLev+1);
                end
            end

            if expShear,
                disCorr = zeros(nF, nR);
                for iLev = 1:nLev-1,
                    inds = pyrBandIndices(pind, iLev);
                    imCorr = phsCmp * (iA2-0.5*(startA+nA-appWrap)) * depPhsRatios{keyLev};
                    disCorr(inds, :) = imresize(imCorr, [length(inds), nR]);
                end
                for iCol = 1:nCol,
                    pyr2(:, :, iCol) = nufftWarp(pyr2(:, :, iCol), -disCorr, pind, filters, 2, 8, keyLev+1, keyLev+1);
                end
            end


            inds = pyrBandIndices(pind, keyLev+1);
            for iCol = 1:nCol,
                errMap(iA - reduWrap, :, :) = squeeze(errMap(iA - reduWrap, :, :)) + imresize(abs(...
                        pyr1(inds, :, iCol) - pyr2(inds, :, iCol)) ./ sqrt(abs(pyr1(inds, :, iCol)).*abs(pyr2(inds, :, iCol))), [nCE, nRE])';
            end
        end
        errMap(isnan(errMap)) = 1e3;
        errMap(errMap > 1e3) = 1e3;
        dispMap = 0;
        nRed = 4;
        for red = nRed:-1:0,
            po2 = 2^red;
            nREr = ceil(nRE/po2);
            nCEr = ceil(nCE/po2);
            cmpMap = zeros(2*(appWrap - reduWrap), nREr, nCEr);
            if isscalar(dispMap),
                for iA = 1:2*(appWrap - reduWrap),
                    cmpMap(iA, :, :) = imresize(squeeze(errMap(iA, :, :)), [nREr, nCEr]);
                end
            else,
                dispMap = imresize(dispMap, [nREr, nCEr]);
                dispMap(dispMap < 1) = 1;
                dispMap(dispMap > 2*(appWrap - reduWrap)) = 2*(appWrap - reduWrap);
                for iA = 1:2*(appWrap - reduWrap),
                    tmp = 1e5*ones(nREr, nCEr);
                    val = imresize(squeeze(errMap(iA, :, :)), [nREr, nCEr]);
                    ind = abs(dispMap - iA) < max(2*(appWrap - reduWrap) / (2^(nRed-red)), 0.6);
                    tmp(ind) = val(ind);
                    cmpMap(iA, :, :) = tmp;
                end
            end
            [unused, dispMap] = min(cmpMap);
            dispMap = squeeze(dispMap);
            mkdir(sprintf('%s/sup/', oprefix));
            imwrite(imresize(0.5 + repmat(dispMap - (appWrap-reduWrap+0.5), [1, 1, 3]) ...
                    / (2*appWrap - 2*reduWrap), [nR, nC]), sprintf('%s/sup/dispMap-%d.jpg', oprefix, red),'jpg');
        end
        dispMap = squeeze(dispMap) - (appWrap-reduWrap+0.5);
    else,
        dispMap = 0;
    end

    for iA = startA:nA - appWrap
        iA - startA + 1
        t0 = cputime();
        img1 = double(imread(ifiles{iA})) / 255;
        pyr1 = zeros([nF, nR, nCol]);

        for iR=1:nR,
            for iCol = 1:nCol,
                pyr1(:, iR, iCol) = buildSCFpyrGen1D(img1(iR, :, iCol), filters, 1, 1, nLev);
            end
        end

        if expShear,
            if false,
                mulCorr = ones(nF, nR);
                for iLev = 1:nLev-1,
                    inds = pyrBandIndices(pind, iLev);
                    imCorr = exp(i * (depPhsRatios{keyLev} / depPhsRatios{iLev}) * phsCmp * (iA-0.5*(startA+nA-appWrap)));
                    mulCorr(inds, :) = imresize(imCorr, [length(inds), nR]);
                end
                mulCorr = mulCorr ./ abs(mulCorr);
                pyr1 = pyr1 ./ repmat(mulCorr, [1, 1, 3]);
            else,
                disCorr = zeros(nF, nR);
                for iLev = 1:nLev-1,
                    inds = pyrBandIndices(pind, iLev);
                    imCorr = phsCmp * (iA-0.5*(startA+nA-appWrap)) * depPhsRatios{keyLev};
                    disCorr(inds, :) = imresize(imCorr, [length(inds), nR]);
                end
                for iCol = 1:nCol,
                    pyr1(:, :, iCol) = nufftWarp(pyr1(:, :, iCol), -disCorr, pind, filters, 2, 8, 1, nLev);
                end
            end
        end

        mixFlag = false;
        if iA + corrDist <= nA,
            mixFlag = true;
            iA2 = iA + corrDist;
        elseif iA - corrDist > 0,
            mixFlag = true;
            iA2 = iA - corrDist;
        end

        if mixFlag,
            img2 = double(imread(ifiles{iA2})) / 255;
            pyr2 = zeros([nF, nR, nCol]);

            for iR=1:nR,
                for iCol = 1:nCol,
                    pyr2(:, iR, iCol) = buildSCFpyrGen1D(img2(iR, :, iCol), filters, 1, 1, nLev);
                end
            end

            if expShear,
                if false,
                    mulCorr = ones(nF, nR);
                    for iLev = 1:nLev-1,
                        inds = pyrBandIndices(pind, iLev);
                        imCorr = exp(i * (depPhsRatios{keyLev} / depPhsRatios{iLev}) * phsCmp * (iA2-0.5*(startA+nA-appWrap)));
                        mulCorr(inds, :) = imresize(imCorr, [length(inds), nR]);
                    end
                    mulCorr = mulCorr ./ abs(mulCorr);
                    pyr2 = pyr2 ./ repmat(mulCorr, [1, 1, 3]);

                else,
                    disCorr = zeros(nF, nR);
                    for iLev = 1:nLev-1,
                        inds = pyrBandIndices(pind, iLev);
                        imCorr = phsCmp * (iA2-0.5*(startA+nA-appWrap)) * depPhsRatios{keyLev};
                        disCorr(inds, :) = imresize(imCorr, [length(inds), nR]);
                    end
                    for iCol = 1:nCol,
                        pyr2(:, :, iCol) = nufftWarp(pyr2(:, :, iCol), -disCorr, pind, filters, 2, 8, 1, nLev);
                    end
                end
            end

            sig = sigWrap * ones(nF, 1);
            dis = zeros(nF, nR);
            for k = nLev:-1:1,
                ind = pyrBandIndices(pind,k);
                sig(ind) = sig(ind) / (2^(nLev - k));%, appWrap/2);
                dis(ind, :) = imresize(dispMap, [nR, length(ind)])';
            end
            sig(sig > appWrap/2) = appWrap/2;
            sig = repmat(sig, [1, nR]);
            if iA - startA < nA-appWrap-iA,
                pref = repmat(0.5 + 0.5 * erf((0.5 + iA-startA - dis)./ (sqrt(2) * sig)), [1, 1, nCol]);
            else,
                pref = repmat(0.5 + 0.5 * erf((0.5 + nA - appWrap-iA + dis)./ (sqrt(2) * sig)), [1, 1, nCol]);
            end
            for k = nLev:-1:1,
                ind = pyrBandIndices(pind,k);
            end
            pyr1 = pyr1 .* pref + pyr2 .* (1 - pref);
        end

        synth = zeros([nR, nC, nCol]);
        for iR = 1:nR,
            for iCol = 1:nCol,
                recons = reconSCFpyrGen1D(pyr1(:, iR, iCol), pind, filters) ;
                synth(iR, :, iCol) = recons;
            end
        end
        imwrite(synth, sprintf('%s/%d.jpg', oprefix, iA-startA), 'jpg');
        cputime() - t0
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

function pyr = buildSCFpyrGen1D(im, filters, expand, minLev, maxLev)
    % Return pyramid in the usual format of a stack of column vectors
    imdft = fftshift(fft(im)); %DFT of image
    nLev = max(size(filters));
    n = length(im);
    pyr = [];
    for k = 1:nLev,
        if k >= minLev && k <= maxLev,
            curFilter = filters{k};
            tempDFT = curFilter.*imdft; % Transform domain

            indices = getIDXFromFilter1D(curFilter, expand);
            nUnex = length(getIDXFromFilter1D(curFilter, 1));
            % shift to avoid cross-boundary wavelet
            tempDFT = tempDFT .* exp(2*pi*1i*((0:n-1) - floor(n/2))/nUnex);
            tempDFT = tempDFT(indices);
            curResult = ifft(ifftshift(tempDFT)) * length(indices) / nUnex;
            curResult = circshift(curResult, [0, floor(length(indices)/(2*nUnex))]);
        else,
            curFilter = filters{k};
            indices = getIDXFromFilter1D(curFilter, expand);
            curResult = zeros(length(indices), 1);
        end
        pyr = [pyr; curResult(:)];
    end
end

function resPyr = nufftWarp(pyr, mov, pind, filters, m, Q, minLev, maxLev)
    %http://www.cscamm.umd.edu/programs/fam04/qing_liu_fam04.pdf
    nLev = max(size(filters));
    nR = size(pyr, 2);
    n = length(filters{1});
    resPyr = zeros(size(pyr));
    for k = maxLev:-1:minLev,
        curFilter = filters{k};
        indices = pyrBandIndices(pind, k);
        if k == nLev,
            mov(indices, :) = 0;
        end
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
        mtFloor = round(mt);
        stack = zeros(nL, 1);
        W = ones(nL, nR);

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
            sig = sig + accumarray(subs, reshape(pyr(indices, :) .* W .* squeeze(X(iQ, :, :)),...
                    [nL*nR, 1]), [m*nL+q, nR]);
            %sig(idx) = sig(idx) + pyr(indices, :) .* W .* squeeze(X(iQ, :, :));
        end
        sig(1+m*nL:q/2+m*nL, :) = sig(1+m*nL:q/2+m*nL, :) + sig(1:q/2, :);
        sig(q/2+1:q, :) = sig(q/2+1:q, :) + sig(q/2+1+m*nL:q+m*nL, :);
        sig = sig(q/2+1:q/2+m*nL, :);

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
