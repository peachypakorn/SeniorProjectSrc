function [synth] = eularian1D(img1, img2, rho, ht, nA, dA, maxDepth, sigWrap)
% clamps depth if maxDepth > 0 using arctan
    nR = size(img1, 1);
    nC = size(img1, 2);
    nCol = size(img1, 3);
    synth = zeros([ size(img1)  nA*2]);
    maxAmp = (nA-1)*dA;

    sizFilt = 2.^[0:-1:-ht];
    filters = get1DRadialFilters(nC, sizFilt, 1);
    nLev = max(size(filters));
    depPhsRatios = getDepthPhaseRatio(filters);

    pind = buildBandIndices(filters);
    nF = pind(nLev, 1) + pind(nLev, 2) - 1;

    pyr1 = zeros([nF nR nCol]);
    pyr2 = zeros([nF nR nCol]);

    for iR=1:nR,
        for iCol = 1:nCol,
            pyr1(:, iR, iCol) = buildSCFpyrGen1D(img1(iR, :, iCol), filters);
            pyr2(:, iR, iCol) = buildSCFpyrGen1D(img2(iR, :, iCol), filters);
        end
    end;
    %% Computing Filter
    d1 = angle(pyr1);
    d2 = angle(pyr2);
    phsDif = (d2 - d1);

    sigma_v_fac = 1;
    phs = exp(phsDif*1i);
    wPhs = zeros([nF, nR]);
    sumPhs = zeros([nF, nR]);
    for c = 1:nCol,
        tmp = bsxfun(@min, abs(pyr1(:, :, c)), abs(pyr2(:, :, c)));
        wPhs = wPhs + tmp;
        sumPhs = sumPhs + phs(:, :, c) .* tmp;
    end
    wPhs(wPhs == 0) = 1e-30;

    for c = 1:nCol,
        for k = nLev:-1:1,
            ind = pyrBandIndices(pind,k);
            cLen = (2^max(0,k-2));
            sigma_v = cLen * sigma_v_fac;
            for d = 1:sigma_v*2,
                weight = bsxfun(@min, abs(pyr1(ind, 1+d:nR, c)), abs(pyr2(ind, 1+d:nR, c))) * ...
                        exp(-(d).^2 / (2 * sigma_v^2));
                sumPhs(ind, 1:nR-d) = sumPhs(ind, 1:nR-d) + phs(ind, 1+d:nR, c) .* weight;
                wPhs(ind, 1:nR-d) = wPhs(ind, 1:nR-d) + weight;

                weight = bsxfun(@min, abs(pyr1(ind, 1:nR-d, c)), abs(pyr2(ind, 1:nR-d, c))) * ...
                        exp(-(d).^2 / (2 * sigma_v^2));
                sumPhs(ind, 1+d:nR) = sumPhs(ind, 1+d:nR) + phs(ind, 1:nR-d, c) .* weight;
                wPhs(ind, 1+d:nR) = wPhs(ind, 1+d:nR) + weight;
            end;
        end
    end

    phsDif = angle(sumPhs ./ wPhs);

    depth = zeros([nF, nR]);
    for k = nLev-1:-1:1,
        ind = pyrBandIndices(pind,k);
        depth(ind, :) = phsDif(ind, :) * depPhsRatios{k};
    end

    phsDifAA = abs(mod(phsDif + 3*pi, 2*pi) - pi);
    for k = nLev-2:-1:1,
        ind = pyrBandIndices(pind,k);
        ind_l = pyrBandIndices(pind,k+1);
        for iR = 1:nR,
            for idx = 1:max(size(ind_l)),
                if k > 1,
                    for d = -1:0,
                        if idx*2+d > 0 && idx*2+d <= max(size(ind)),
                            if phsDifAA(ind_l(idx), iR) >= pi / 2,
                                %phsDifAA(ind(idx*2+d), iR) <= phsDifAA(ind_l(idx), iR) * 2 && ...
                                %    phsDifAA(ind_l(idx), iR) * (1+maxAmp) >= pi / 2,
                                phsDifAA(ind(idx*2+d), iR) = phsDifAA(ind_l(idx), iR) * 2;
                                depth(ind(idx*2+d), iR) = depth(ind_l(idx), iR);
                            end
                        end
                    end
                else,
                    if phsDifAA(ind_l(idx), iR) >= pi / 2,
                        %phsDifAA(ind(idx*2+d), iR) <= phsDifAA(ind_l(idx), iR) * 2 && ...
                        %    phsDifAA(ind_l(idx), iR) * (1+maxAmp) >= pi / 2,
                        phsDifAA(ind(idx), iR) = phsDifAA(ind_l(idx), iR) * 2;
                        depth(ind(idx), iR) = depth(ind_l(idx), iR);
                    end
                end
            end
        end
    end

    if maxDepth > 0,
        depthClp = maxDepth * atan(depth / maxDepth);
        %depthClp = clip(depth, -maxDepth, maxDepth);
    else
        depthClp = depth;
    end
    %for k = nLev-1:-1:1,
    %    ind = pyrBandIndices(pind,k);
    %    indChg = abs(depthClp(ind) / depPhsRatios{k}) > pi/(2*nA*dA);
    %    depthClp(ind(indChg)) = round(depthClp(ind(indChg)) * 2*nA*dA / ...
    %            (depPhsRatios{k} * 2 * pi)) * (pi/(nA*dA)) * depPhsRatios{k};
    %end

    for k = nLev-1:-1:1,
        ind = pyrBandIndices(pind,k);
        phsDifAA(ind, :) = depthClp(ind, :) / depPhsRatios{k};
    end

    for c = 1:size(img1,3),
        %% Reconstructing
        appWrap = ceil(sigWrap*2);
        pyrRec = zeros(nF, nR, 2*nA+2*appWrap);
        for iA = 1:nA+appWrap,
            amp = -0.5 + (iA - 0.5)*dA;

            depChg = -depthClp * (0.5+amp) + depth/2;
            depChg2 = -depthClp * (0.5+amp) - depth/2;
            phsChg = zeros(size(depChg));
            phsChg2 = zeros(size(depChg2));
            fcc = abs(depChg);
            fcc(fcc == 0) = abs(amp*1e-10);
            fcc2 = abs(depChg2);
            fcc2(fcc2 == 0) = abs((1+amp)*1e-10);
            mixFactor = fcc2 ./ (fcc + fcc2);
            for k = nLev-1:-1:1,
                ind = pyrBandIndices(pind,k);
                phsChg(ind, :) = depChg(ind, :) / depPhsRatios{k};
                phsChg2(ind, :) = depChg2(ind, :) / depPhsRatios{k};
            end
            if iA > nA - appWrap,
                fRho = 1 + (iA - nA + appWrap) * nA / appWrap;
            else,
                fRho = 1;
            end

            pyrRec(:, :, nA+appWrap+1-iA) = exp( -( phsDifAA * rho *fRho).^2 / 2).* (...
                    nufftWarp(pyr1(:, :, c), -depChg, -depthClp, pind, filters, 2, 8) .* mixFactor +...
                    nufftWarp(pyr2(:, :, c), -depChg2, -depthClp, pind, filters, 2, 8) .* (1-mixFactor));

            pyrRec(:, :, nA+appWrap+iA) = exp( -( phsDifAA * rho *fRho).^2 / 2).* (...
                    nufftWarp(pyr2(:, :, c), depChg, -depthClp, pind, filters, 2, 8) .* mixFactor +...
                    nufftWarp(pyr1(:, :, c), depChg2, -depthClp, pind, filters, 2, 8) .* (1-mixFactor));
            %pyrRec(:, :, nA+appWrap+1-iA) = exp( -( phsDifAA * rho *fRho).^2 / 2).* (...
            %        exp( 1i* phsChg) .* pyr1(:, :, c) .* mixFactor +...
            %        exp( 1i* phsChg2) .* pyr2(:, :, c) .* (1-mixFactor));

            %pyrRec(:, :, nA+appWrap+iA) = exp( -( phsDifAA * rho *fRho).^2 / 2).* (...
            %        exp( 1i* -phsChg) .* pyr2(:, :, c) .* mixFactor +...
            %        exp( 1i* -phsChg2) .* pyr1(:, :, c) .* (1-mixFactor));
        end

        for iA = 1:appWrap,
            mixFactor = exp((iA - 0.5) * sigWrap) / (1 + exp((iA - 0.5) * sigWrap));
            pyrRec(:, :, appWrap + iA) = pyrRec(:, :, appWrap + iA) * mixFactor + ...
                    pyrRec(:, :, 2*nA + appWrap + iA) * (1-mixFactor);
            pyrRec(:, :, 2*nA + appWrap + 1 - iA) = pyrRec(:, :, 2*nA + appWrap + 1 - iA) * mixFactor + ...
                    pyrRec(:, :, appWrap + 1 - iA) * (1-mixFactor);
        end

        for iA = 1:2*nA,
            for iR = 1:nR,
                recons = reconSCFpyrGen1D(pyrRec(:, iR, iA+appWrap), pind, filters) ;
                synth(iR, :, c, iA) = recons;
            end
        end
    end
end

function outFilters = get1DRadialFilters(dims, rVals, twidth)
    % Construct log_rad
    ctr = ceil((dims+0.5)/2);
    log_rad = abs(((1:dims)-ctr)./(dims/2));
    log_rad(ctr) = log_rad(ctr-1);

    N = max(size(rVals));

    [himask, lomaskPrev] = getRadialMaskPair(rVals(1),log_rad, twidth);
    outFilters{1} = himask;
    for k = 2:N
        [himask, lomask] = getRadialMaskPair(rVals(k),log_rad, twidth);
        mask = (1:dims)>ctr;
        outFilters{k} = himask.*lomaskPrev.*mask;
        lomaskPrev = lomask;
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

function pind = buildBandIndices(filters)
    nLev = max(size(filters));
    pind = [];
    sumIdx = 1;
    for k = 1:nLev
        curFilter = filters{k};
        indices = getIDXFromFilter1D(curFilter);
        pind = [pind; sumIdx, numel(indices)];
        sumIdx = sumIdx + numel(indices);
    end
end

function indices =  pyrBandIndices(pind,band)
    indices = pind(band, 1):pind(band, 1) + pind(band, 2) - 1;
end

function pyr = buildSCFpyrGen1D(im, filters)
    % Return pyramid in the usual format of a stack of column vectors
    imdft = fftshift(fft(im)); %DFT of image
    nLev = max(size(filters));
    n = length(im);
    pyr = [];
    for k = 1:nLev
        curFilter = filters{k};
        tempDFT = curFilter.*imdft; % Transform domain

        indices = getIDXFromFilter1D(curFilter);
        % shift to avoid cross-boundary wavelet
        tempDFT = tempDFT .* exp(2*pi*1i*((0:n-1) - floor(n/2))/length(indices));

        tempDFT = tempDFT(indices);
        curResult = ifft(ifftshift(tempDFT));
        pyr = [pyr; curResult(:)];
    end
end

function resPyr = nufftWarp(pyr, mov, depth, pind, filters, m, Q)
    %http://www.cscamm.umd.edu/programs/fam04/qing_liu_fam04.pdf
    nLev = max(size(filters));
    nR = size(pyr, 2);
    n = length(filters{1});
    resPyr = zeros(size(pyr));
    for k = 1:nLev,
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
        mtFloor = round(mt);
        stack = zeros(nL);
        for iR = 1:nR,
            top = 1;
            for iL = 1:nL,
                while top > 1 && mtFloor(iL, iR) <= mtFloor(stack(top-1), iR) ...
                        && depth(indices(iL), iR) > depth(indices(stack(top-1)), iR),
                    pyr(indices(stack(top-1)), iR) = 0;
                    top = top - 1;
                end
                if top == 1 || mtFloor(iL, iR) > mtFloor(stack(top-1), iR),
                    stack(top) = iL;
                    top = top + 1;
                else,
                    pyr(indices(iL), iR) = 0;
                end
            end
        end
        mtFloor = mod(mtFloor, nL*m);

        for iQ = -q/2:q/2,
            frac = mt - mtFloor - iQ;
            tmp = -1i*(sin((pi/m)*(frac-0.5))./(1-exp(1i*2*pi*(frac-0.5)/(nL*m))) +...
                       sin((pi/m)*(frac+0.5))./(1-exp(1i*2*pi*(frac+0.5)/(nL*m))));
            A(iQ+q/2+1, :, :) = tmp;
        end
        X = reshape(F * reshape(A, [q+1, nL*nR]), [q+1, nL, nR]);
        sig = zeros(m*nL+q, nR);
        for iQ = 1:q+1,
            %assumes uniqueness in mtFloor
            idx = mtFloor+iQ + repmat((0:nR-1)*(m*nL+q), [nL, 1]);
            sig(idx) = sig(idx) + pyr(indices, :) .* squeeze(X(iQ, :, :));
        end
        sig(1+m*nL:q/2+m*nL, :) = sig(1+m*nL:q/2+m*nL, :) + sig(1:q/2, :);
        sig(q/2+1:q, :) = sig(q/2+1:q, :) + sig(q/2+1+m*nL:q+m*nL, :);
        sig = sig(q/2+1:q/2+m*nL, :);

        for iR = 1:nR,
            fftSig = fftshift(fft(sig(:, iR)));
            fftSig = fftSig .* cos(1i*pi*((0:m*nL-1)-m*nL/2)/(m*nL))';
            resPyr(indices, iR) = ifft(ifftshift(fftSig(ceil(m*nL/4)+(1:m*nL/2))));
        end
    end
end

function res = reconSCFpyrGen1D(pyr, pind, filters)
    imdft = filters{1}.*fftshift(fft(pyr(pyrBandIndices(pind,1))')); %DFT of image
    nLev = max(size(filters));
    n = length(filters{1});
    for k = 2:nLev
        curFilter = filters{k};
        tempDFT = fftshift(fft(2*real(pyr(pyrBandIndices(pind,k))'))); % Transform domain
        if k < nLev,
            curFilter = curFilter + rot90(curFilter,2);
        else
            tempDFT = tempDFT / 2;
        end
        indices = getIDXFromFilter1D(curFilter);
        imdft(indices) = imdft(indices) + tempDFT.*curFilter(indices)...
                .*exp(-2*pi*1i*(indices-1-floor(n/2))/length(indices));
    end
    res =  real(ifft(ifftshift(imdft)));
end

function [ out ] = getIDXFromFilter1D(filt)
    aboveZero = filt>1e-10;
    aboveZero = logical(aboveZero + fliplr(aboveZero));

    idx1 = 1:length(filt);
    idx1 = idx1(aboveZero);
    idx1 = min(idx1):max(idx1);

    out = idx1;
end