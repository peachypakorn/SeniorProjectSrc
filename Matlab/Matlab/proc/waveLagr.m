function [synth, depthLo, depthRo] = waveLagr(img1, img2, rho, ht, nA, dA, dsrLagr, maxDLagr, dReg, maxDepth, sigWrap, recR, depthLo, depthRo)
% clamps depth if maxDepth > 0 using arctan
    nR = size(img1, 1);
    nC = size(img1, 2);
    nCol = size(img1, 3);
    synth = zeros([ size(img1)  nA*2]);
    maxAmp = (nA-1)*dA;
    sigma_v_fac = 1;

    if isscalar(depthLo) || isscalar(depthRo),
        [depthLo, depthRo] = lagrDepthEstimate(imresize(img1, dsrLagr), imresize(img2, dsrLagr), maxDLagr, dReg);
    end
    lagrDepthL = depthLo;
    lagrDepthR = depthRo;

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

    depth1 = zeros([nF nR]);
    depth2 = zeros([nF nR]);

    for iLev = 1:nLev,
        ind = pyrBandIndices(pind, iLev);
        nO = pind(iLev, 2);

        depth1(ind, :) = imresize(lagrDepthL, [nR, pind(iLev, 2)])' / dsrLagr;
        depth2(ind, :) = imresize(lagrDepthR, [nR, pind(iLev, 2)])' / dsrLagr;
    end

    for iR=1:nR,
        for iCol = 1:nCol,
            pyr1(:, iR, iCol) = buildSCFpyrGen1D(img1(iR, :, iCol), filters, 1);
            pyr2(:, iR, iCol) = buildSCFpyrGen1D(img2(iR, :, iCol), filters, 1);
        end
    end;

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
    end
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
            fcc11(fcc11 == 0) = abs(amp*1e-10);
            fcc12 = abs(depChg12);
            fcc12(fcc12 == 0) = abs((1+amp)*1e-10);
            fcc21 = abs(depChg21);
            fcc21(fcc21 == 0) = abs(amp*1e-10);
            fcc22 = abs(depChg22);
            fcc22(fcc22 == 0) = abs((1+amp)*1e-10);
            %mixFactor1 = fcc22 ./ (fcc11 + fcc22);
            %mixFactor2 = fcc12 ./ (fcc12 + fcc21);
            mixFactor1 = ones(size(fcc11));
            mixFactor2 = ones(size(fcc11));
            for k = nLev:-1:1,
                ind = pyrBandIndices(pind,k);
                phsChg11(ind, :) = depChg11(ind, :) / depPhsRatios{k};
                phsChg12(ind, :) = depChg12(ind, :) / depPhsRatios{k};
                phsChg21(ind, :) = depChg21(ind, :) / depPhsRatios{k};
                phsChg22(ind, :) = depChg22(ind, :) / depPhsRatios{k};
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
end

function [depthL, depthR] = lagrDepthEstimate(Il, Ir, maxDLagr, dReg)
    % https://www.ims.tuwien.ac.at/publications/tuw-210567

    % Parameter settings
    r = 9;                  % filter kernel in eq. (3) has size r \times r
    eps = 0.0001;           % \epsilon in eq. (3)
    thresColor = 7/255;     % \tau_1 in eq. (5)
    thresGrad = 2/255;      % \tau_2 in eq. (5)
    gamma = 0.11;           % (1- \alpha) in eq. (5)
    threshBorder = 3/255;   % some threshold for border pixels
    gamma_c = 0.1;          % \sigma_c in eq. (6)
    gamma_d = 9;            % \sigma_s in eq. (6)
    r_median = 19;          % filter kernel of weighted median in eq. (6) has size r_median \times r_median

    % Convert to grayscale
    Il_g = rgb2gray(Il);
    Ir_g = rgb2gray(Ir);

    % Mirror image
    Il_1 = flipdim(Ir,2);
    Ir_1 = flipdim(Il,2);

    % Compute gradient in X-direction from grayscale images
    fx_l = gradient(Il_g);
    fx_r = gradient(Ir_g);
    fx_l = fx_l+0.5; % To get a range of values between 0 to 1
    fx_r = fx_r+0.5; % To get a range of values between 0 to 1
    fx_l_1 = flipdim(fx_r,2);
    fx_r_1 = flipdim(fx_l,2);

    [m,n,c] = size(Il);

    dispVol = ones(m,n,2*maxDLagr+1)*(threshBorder);
    dispVol1 = ones(m,n,2*maxDLagr+1)*(threshBorder);

    % Create initial cost volume (eq. (5) in the paper)
    for d=-maxDLagr:maxDLagr
        if d > 0,
            range1 = d+1:n;
            range2 = 1:n-d;
        else,
            range1 = 1:n+d;
            range2 = -d+1:n;
        end

        % Truncated SAD of color images for current displacement
        tmp = ones(m,n,c)*threshBorder;
        tmp(:,range1,:) = Ir(:,range2,:);
        p_color = abs(tmp - Il);
        p_color = sum(p_color,3)*0.333333333333;
        p_color = min(p_color,thresColor);

        % Truncated SAD of gradient images for current displacement
        tmp = ones(m,n)*threshBorder;
        tmp(:,range1) = fx_r(:,range2);
        p_grad = abs(tmp - fx_l);
        p_grad = min(p_grad,thresGrad);

        p = gamma*p_color + (1-gamma)*p_grad; % Combined color and gradient

        % Same for other view
        tmp1 = ones(m,n,c)*threshBorder;
        tmp1(:,range1,:) = Ir_1(:,range2,:);
        p1_color = abs(tmp1 - Il_1);
        p1_color = sum(p1_color,3)*0.333333333333;
        p1_color = min(p1_color,thresColor);

        tmp1 = ones(m,n)*threshBorder;
        tmp1(:,range1) = fx_r_1(:,range2);
        p1_grad = abs(tmp1 - fx_l_1);
        p1_grad = min(p1_grad,thresGrad);

        p1 = gamma*p1_color + (1-gamma)*p1_grad; % Combined color and gradient

        % Set value in cost volume
        dispVol(:,:,d+1+maxDLagr) = p;
        dispVol1(:,:,d+1+maxDLagr) = flipdim(p1,2);
    end

    % Smooth cost volume with guided filter (using eq. (1) with weights in (4))
    for d=1:2*maxDLagr+1 % use regular for loop when not using the parallel computing toolbox
        p = dispVol(:,:,d);
        p1 = dispVol1(:,:,d);

        q = guidedfilter_color(Il, double(p), r, eps);
        p1 =  flipdim(p1,2);
        q1 = guidedfilter_color(Il_1, double(p1), r, eps);

        dispVol(:,:,d) = q * (1 + ((d - maxDLagr - 1) / dReg) ^2);
        dispVol1(:,:,d) = flipdim(q1,2) * (1 + ((d - maxDLagr - 1) / dReg) ^2);
    end

    % Winner take all label selection (eq. (2))
    [unused,labels_left] = min(dispVol,[],3);
    [unused,labels_right] = min(dispVol1,[],3);
    labels_left = labels_left - (1+maxDLagr);
    labels_right = labels_right - (1+maxDLagr);

    %figure
    %ana = zeros([size(labels_left), 3]);
    %ana(:, :, 1) = 0.5*labels_left / maxDLagr + 0.5;
    %ana(:, :, 2:3) = repmat(0.5*labels_right / maxDLagr + 0.5, [1, 1, 2]);
    %imshow(ana);

    % Left-right consistency check
    Y = repmat((1:m)', [1 n]);
    X = repmat(1:n, [m 1]) - labels_left;
    X(X<1) = 1;
    X(X>n) = n;
    indices = sub2ind([m,n],Y,X);

    final_labels_left = labels_left;
    final_labels_left(abs(labels_left - labels_right(indices))>=1) = -maxDLagr-1;

    %figure
    %ana(:, :, 1) = 0.5*final_labels_left / maxDLagr + 0.5;

    % Fill and filter (post-process) pixels that fail the consistency check
    final_labels_left = fillPixelsReference(Il, final_labels_left, gamma_c, gamma_d, r_median, maxDLagr);

    depthL = final_labels_left;

    % Left-right consistency check
    Y = repmat((1:m)', [1 n]);
    X = repmat(1:n, [m 1]) + labels_right;
    X(X<1) = 1;
    X(X>n) = n;
    indices = sub2ind([m,n],Y,X);

    final_labels_right = labels_right;
    final_labels_right(abs(labels_right - labels_left(indices))>=1) = -maxDLagr-1;

    %ana(:, :, 2:3) = repmat(0.5*final_labels_right / maxDLagr + 0.5, [1, 1, 2]);
    %imshow(ana);

    % Fill and filter (post-process) pixels that fail the consistency check
    final_labels_right = fillPixelsReference(Ir, final_labels_right, gamma_c, gamma_d, r_median, maxDLagr);

    depthR = final_labels_right;
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
            sig = sig + accumarray(subs, reshape(pyr(indices, :) .* W .* squeeze(X(iQ, :, :)),...
                    [nL*nR, 1]), [m*nL+q, nR]);
            %sig(idx) = sig(idx) + pyr(indices, :) .* W .* squeeze(X(iQ, :, :));
            if k == nLev,
                sig2 = sig2 + accumarray(subs, reshape(W .* squeeze(X(iQ, :, :)),...
                        [nL*nR, 1]), [m*nL+q, nR]);

                %sig2(idx) = sig2(idx) + W .* squeeze(X(iQ, :, :));
            end
        end
        sig(1+m*nL:q/2+m*nL, :) = sig(1+m*nL:q/2+m*nL, :) + sig(1:q/2, :);
        sig(q/2+1:q, :) = sig(q/2+1:q, :) + sig(q/2+1+m*nL:q+m*nL, :);
        sig = sig(q/2+1:q/2+m*nL, :);
        if k == nLev,
            sig2(1+m*nL:q/2+m*nL, :) = sig2(1+m*nL:q/2+m*nL, :) + sig2(1:q/2, :);
            sig2(q/2+1:q, :) = sig2(q/2+1:q, :) + sig2(q/2+1+m*nL:q+m*nL, :);
            sig2 = sig2(q/2+1:q/2+m*nL, :);
        end

        for iR = 1:nR,
            fftSig = fftshift(fft(sig(:, iR)));
            fftSig = fftSig .* cos(1i*pi*((0:m*nL-1)-m*nL/2)/(m*nL))';
            resPyr(indices, iR) = ifft(ifftshift(fftSig(ceil(m*nL/4)+(1:m*nL/2))));
            if k == nLev,
                fftSig2 = fftshift(fft(sig2(:, iR)));
                fftSig2 = fftSig2 .* cos(1i*pi*((0:m*nL-1)-m*nL/2)/(m*nL))';
                dem = ifft(ifftshift(fftSig2(ceil(m*nL/4)+(1:m*nL/2))));
                dem(dem == 0) = 1;
                resPyr(indices, iR) = resPyr(indices, iR) ./ dem;
            end
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
