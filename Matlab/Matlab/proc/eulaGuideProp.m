function [synth, mem] = eulaOnLagrGuideProp(img1, img2, testParam, runParam, mem, iA, phsCmpO)
    rho = runParam.rho;
    ht = runParam.ht;
    nA = runParam.nA;
    dA = runParam.dA;
    
    if isfield(runParam, 'remap'),
        maxDepth = runParam.remapDepth;
    else
        maxDepth = 0;
    end
    if isfield(runParam, 'recR'),
        recR = runParam.recR;
    else,
        recR = 0;
    end

    waveletExt = 4;
    nR = size(img1, 1);
    nC = size(img1, 2);
    nCol = size(img1, 3);

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
        rec1 = zeros([nF nR nCol]);
        rec2 = zeros([nF nR nCol]);
        rec1Ref = zeros([nFRef nR nCol]);
        rec2Ref = zeros([nFRef nR nCol]);

        refInd1 = ones([nF nR nCol]);
        refInd2 = ones([nF nR nCol]);
        preDepth1 = zeros([nF nR]);
        preDepth2 = zeros([nF nR]);

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
        tempDFT = tempDFT(indices);
        curResult = ifft(ifftshift(tempDFT)) * length(indices) / nUnex;
        pyr = [pyr; curResult(:)];
    end
end

function pyr = buildSCFpyrGen1DSum(im, filters, expand)
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
        tempDFT = tempDFT(indices);
        curResult = ifft(ifftshift(tempDFT)) * length(indices) / nUnex;
        pyr = [pyr; curResult(:)];
    end
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