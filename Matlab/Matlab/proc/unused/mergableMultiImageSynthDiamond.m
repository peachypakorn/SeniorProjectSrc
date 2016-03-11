function [synth, depthE, movPyr] = mergableMultiPyrSynth(img, imgref, amp, rStripDebug, trim, dscale, doffset, rho, ...
        rhor, mmtrim, prvMovPyr, prvAlpha)
    cmpSig = zeros([size(img, 1), size(img, 2)]);
    cmpRef = zeros([size(img, 1), size(img, 2)]);

    for r = 1:size(img, 1),
        bw = rgb2ntsc(squeeze(img(r, :, :)));
        cmpSig(r, :) = bw(:, 1);
        bw = rgb2ntsc(squeeze(imgref(r, :, :)));
        cmpRef(r, :) = bw(:, 1);
    end
    synth = zeros([ size(img) length(amp) ]);
    depthP = zeros([ size(img) 1+ceil(log2(size(img, 2) / trim))]);
    depth = zeros([ size(img) ]);
    depthE = zeros([ size(img) ]);

    t0 = cputime();
    [movPyr, dep] = movePyrRec(cmpSig, cmpRef, trim);
    cputime() - t0

    depthE = repmat(dep, [1, 1, 3]);
    if rStripDebug > 0,
        out_folder = sprintf('result/tmp_debug');
        mkdir(out_folder);
        imwrite(depthE(:, :, :) / dscale + doffset, sprintf('%s/depthE.png', out_folder),'png');
    end

    if prvAlpha ~= 0 & iscell(prvMovPyr) & length(movPyr) == length(prvMovPyr),
        for lay = 1:length(movPyr),
            gam = prvAlpha * size(img, 2) / (2^(lay-1));
            w = exp(-((movPyr{lay} - prvMovPyr{lay})/gam).^2);
            movPyr{lay} = movPyr{lay} .* (1 ./ (1+w)) + prvMovPyr{lay} .* (w ./ (1+w));
        end
    end

    'pyrimid'  
    t0 = cputime();
    for r = 1:size(img, 1),
        for c = 1:3,
            synth(r, :, c, :) = mergableMultiPyrFast(img(r, :, c), cmpSig(r, :), cmpRef(r, :), movPyr, r, trim, rho, rhor, amp, mmtrim);
            %[synth(r, :, c, :), dep] = mergableMultiPyrRec(img(r, :, c), movPyr, r, trim, rho,...
            %        amp, 1, 1, 1, size(img, 2), ones([1 size(img, 2)]), 0);
            %for l = 1:length(dep),
            %    if l < length(dep),
            %       depthP(r, :, :, length(dep)+1-l) = repmat(dep{l}' - dep{l+1}', 1, 3);
            %    else,
            %       depthP(r, :, :, length(dep)+1-l) = repmat(dep{l}', 1, 3);
            %    end
            %end
            %depth(r, :, :) = repmat(dep{1}', 1, 3);
        end
        %r
        if rStripDebug > 0,
            if mod(r, rStripDebug) == 0,
                r
                %for l = 1:length(dep),
                %    imwrite(depthP(:, :, :, l) / min(dscale, size(img, 1) / 2^(l-1)) + 0.5, sprintf('%s/depth%03d.png', out_folder, l),'png');
                %end
                %imwrite(depth(:, :, :) / dscale + doffset, sprintf('%s/depth.png', out_folder),'png');
                for m=1:length(amp),
                    imwrite(repmat(synth(:, :, :, m), [1, 1, 1]), sprintf('%s/%03d.png', out_folder, m+2),'png');
                end
            end
        end
    end
    cputime() - t0
end

function [movPyr, depth] = movePyrRec(cmpSig, cmpRef, trim)
    nC = size(cmpSig, 2);
    nR = size(cmpSig, 1);
    powW = 2.0;
    for cLay = 1:round(log2(nC / trim)) + 1,
        movPyr{cLay} = zeros(nR, 2^cLay - 1);
        cLen = nC * 2 / (2^cLay)
        filt = 1:cLen;
        filt = cos((filt - (cLen+1)/2) * pi / cLen);

        samp = min(cLen/4, 16);
        xx = 1:samp;
        moves = zeros(1, 8);
        freqSig = fft(cmpSig, [], 2) ;
        freqRef = fft(cmpRef, [], 2) ;
        filtSig = zeros([nR*nC samp]);
        filtRef = zeros([nR*nC samp]);
        for s = 1:samp,
            freqFilt = zeros(1, nC);
            freqFilt(nC:-1:nC-cLen+1) = filt .* exp(-2 * pi * 1i * (0:cLen-1) * s / cLen);
            freqFilt(nC:-1:nC-cLen+1) = freqFilt(nC:-1:nC-cLen+1) - filt * sum(freqFilt(nC:-1:nC-cLen+1)) / sum(filt);
            freqFilt = fft(freqFilt);
            filtSig(:, s) = reshape(ifft(freqSig .* repmat(freqFilt, nR, 1), [], 2), [nR*nC, 1]);
            filtRef(:, s) = reshape(ifft(freqRef .* repmat(freqFilt, nR, 1), [], 2), [nR*nC, 1]);
        end
        for r = 1:nR,
            r1 = clip(r - cLen / 2, 1, nR);
            r2 = clip(r + cLen / 2, 1, nR);
            for c = 1:2^cLay - 1,
                nM = 0;
                if cLay == 1,
                    moves(1) = 0;
                    nM = 1;
                else,
                    p = floor(c / 2 + 1e-10);
                    if p > 0 && p <= 2^(cLay - 1) - 1,
                        nM = nM + 1;
                        moves(nM) = movPyr{cLay - 1}(r, p);
                        nM = nM + 1;
                        moves(nM) = movPyr{cLay - 1}(r1, p);
                        nM = nM + 1;
                        moves(nM) = movPyr{cLay - 1}(r2, p);
                    end
                    p = floor((c + 1) / 2 + 1e-10);
                    if p > 0 && p <= 2^(cLay - 1) - 1 && p ~= floor(c / 2 + 1e-10),
                        nM = nM + 1;
                        moves(nM) = movPyr{cLay - 1}(r, p);
                        nM = nM + 1;
                        moves(nM) = movPyr{cLay - 1}(r1, p);
                        nM = nM + 1;
                        moves(nM) = movPyr{cLay - 1}(r2, p);
                    end
                    p = floor(c / 2 + 1e-10) - 1;
                    if p > 0 && p <= 2^(cLay - 1) - 1,
                        nM = nM + 1;
                        moves(nM) = movPyr{cLay - 1}(r, p);
                    end
                    p = floor((c + 1) / 2 + 1e-10) + 1;
                    if p > 0 && p <= 2^(cLay - 1) - 1,
                        nM = nM + 1;
                        moves(nM) = movPyr{cLay - 1}(r, p);
                    end
                end
                st = (c - 1) * cLen / 2 + 1;
                sig = filtSig(r + (st - 1) * nR, :);
                moves = round(moves + st);
                if all(moves <= -cLen / 8) || all(moves > nC - cLen + 1 + cLen / 8),
                    movPyr{cLay}(r, c) = moves(1) - st;
                else
                    moves = clip(moves, 1, nC - cLen + 1);
                    sumW = 0;
                    sumMov = 0;
                    for mId = 1:nM,
                        rst = moves(mId);
                        ref = filtRef(r + (rst - 1) * nR, :);
                        ang = angle(sig ./ (ref + 1e-100));
                        weight = abs(sig .* ref);
                        wx = weight .* xx;
                        wx2 = weight .* (xx .^2);

                        id = 4;
                        while id <= samp,
                            angM = sum(ang(1:id).*wx(1:id)) / (sum(wx2(1:id)) + 1e-100);
                            ang = ang + round((angM * (1:samp) - ang) / (2*pi)) * 2 * pi;
                            id = id + 4;
                        end
                        angM = sum(ang .* wx) / (sum(wx2) + 1e-100);
                        err = norm((angM * xx - ang) .* weight) / (sum(weight) + 1e-100) + 1e-30;
                        nMov = rst - st + angM * cLen / (2*pi);
                        sumMov = sumMov + nMov / err ^ powW;
                        sumW = sumW + 1.0 / err^powW;
                    end
                    movPyr{cLay}(r, c) = sumMov / sumW;
                end
            end
        end

        filtered = zeros(size(movPyr{cLay}));
        rSigma = cLen / 3;
        hKernel = ceil(rSigma * 1.5);
        dSigma = cLen / 2;
        for r = 1:nR,
            for c = 1:2^cLay - 1,
                w = 0;
                for dr = -hKernel:hKernel,
                    if dr == 0,
                        continue;
                    end
                    r2 = r + dr;
                    while r2 <= 0 || r2 > nR,
                        if r2 <= 0,
                            r2 = 1 - r2;
                        elseif r2 > nR,
                            r2 = 2*nR + 1 - r2;
                        end
                    end
                    dw = exp(-((r2 - r)^2/rSigma^2 + (movPyr{cLay}(r2, c) - movPyr{cLay}(r, c))^2/dSigma^2));
                    filtered(r, c) = filtered(r, c) + dw * movPyr{cLay}(r2, c);
                    w = w + dw;
                end
                filtered(r, c) = filtered(r, c) / w;
            end
        end
        movPyr{cLay} = filtered;
    end
    depth = zeros(nR, nC);
    cLay = round(log2(nC / trim)) + 1;
    cLen = nC * 2 / (2^cLay);
    filt = 1:cLen;
    filt = cos((filt - (cLen+1)/2) * pi / cLen);
    for r = 1:nR,
        for c = 1:2^cLay - 1,
            st = (c - 1) * cLen / 2 + 1;
            depth(r, st:st+cLen-1) = depth(r, st:st+cLen-1) + (filt.^2) * movPyr{cLay}(r, c);
        end
    end
end

function synth = mergableMultiPyrFast(sig, cmpSig, cmpRef, movPyr, r, trim, rho, rhor, amp, mmtrim)
    nC = length(sig);
    nA = length(amp);
    nLay = round(log2(nC / trim)) + 1;
    synth = zeros(nC, nA);
    for cLay = 1:nLay,
        cLen = nC * 2 / (2^cLay);
        sPyr{cLay} = zeros(2^cLay - 1, cLen, nA);
    end

    mmfilters = [];
    if mmtrim < trim & mmtrim ~= 0,
        mmLayers = round(log2(trim / mmtrim));
        mmfilters = zeros(mmLayers+1, trim);
        lg = zeros(1, trim);
        lg(2:trim/2) = log2((1:trim/2-1) / (trim/2));

        res = ones(1, trim);
        for mLay = 1:mmLayers,
            hp = -sin(clip((lg - log2(2/trim)) - mLay, -1, 0) * pi / 2);
            mmfilters(mLay, :) = res .* hp;
            res = res .* sqrt(1 - hp.^2);
        end
        mmfilters(mmLayers+1, :) = res;
    end

    for cLay = nLay:-1:1,
        cLen = nC * 2 / (2^cLay);
        filt = 1:cLen;
        filt = cos((filt - (cLen+1)/2) * pi / cLen);
        lFilt = 1:cLen*2;
        lFilt = cos((lFilt - (2*cLen+1)/2) * pi / (2*cLen));
        
        dif = 0:cLen-1;
        dif(cLen/2+1) = 0;
        dif(cLen:-1:cLen/2+2) = -1:-1:-cLen/2+1;

        if cLay ~= nLay && cLay ~= 1 % only an estimation
            fac = abs(dif);
            fac(4:2:cLen-2) = 2*fac(4:2:cLen-2);
            fac2 = dif .^ 2 - fac;
        elseif cLay == nLay,
            fac = abs(dif);
            fac(4:2:cLen-2) = 2*fac(4:2:cLen-2);
            fac2 = 0;
        elseif cLay == 1,
            fac = abs(dif);
            fac(4:2:cLen-2) = 2 * fac(4:2:cLen-2);
            fac2 = 0;
        end

        for c = 1:2^cLay-1,
            cSt = (c - 1) * cLen / 2 + 1;
            iFilt1 = 1;
            iFilt2 = 1;
            if cLay == 1,
                cFilt = ones(size(filt));
            elseif c == 1,
                cFilt = ones(size(filt));
                cFilt(cLen/2+1:cLen) = filt(cLen/2+1:cLen);
            elseif c == 2^cLay-1,
                cFilt = ones(size(filt));
                cFilt(1:cLen/2) = filt(1:cLen/2);
            else,
                cFilt = filt;
                iFilt1 = ones(size(filt));
                iFilt1(1:cLen/2) = filt(1:cLen/2);
                iFilt2 = ones(size(filt));
                iFilt2(cLen/2+1:cLen) = filt(cLen/2+1:cLen);
            end

            if cLay == nLay,
                cSt = (c - 1) * cLen / 2 + 1;
                if mmtrim >= trim | mmtrim == 0,
                    sPyr{cLay}(c, :, :) = repmat(sig(cSt:cSt+cLen-1) .* cFilt, [1, 1, nA]);
                else
                    sPyr{cLay}(c, :, :) = motionAmp(sig, cmpSig, cmpRef, cSt, trim, movPyr{cLay}(r, c), amp, rho, mmfilters, filt);
                end
            end
            pars = [floor(c / 2 + 1e-10), floor((c + 1) / 2 + 1e-10)];
            if pars(1) == pars(2),
                pars = [pars(1)];
            end

            for pId = 1:length(pars),
                p = pars(pId);
                if cLay == 1,
                    mov = 0;
                elseif p <= 0 || p > 2^(cLay - 1) - 1,
                    continue;
                else,
                    mov = movPyr{cLay - 1}(r, p);
                end
                if cLay == 2,
                    pFilt = ones(size(lFilt));
                elseif p == 1,
                    pFilt = ones(size(lFilt));
                    pFilt(cLen+1:2*cLen) = lFilt(cLen+1:2*cLen);
                elseif p == 2^(cLay-1)-1,
                    pFilt = ones(size(lFilt));
                    pFilt(1:cLen) = lFilt(1:cLen);
                else,
                    pFilt = lFilt;
                end

                nMov = movPyr{cLay}(r, c) - mov;
                for id = 1:nA,
                    if p * 2 - 1 == c,
                        cSig = sPyr{cLay}(c, :, id) .* pFilt(1:cLen) ./ iFilt1;
                        aFilt = cFilt .* pFilt(1:cLen) ./ iFilt1;
                    elseif p * 2 == c,
                        cSig = sPyr{cLay}(c, :, id) .* pFilt(cLen/2+1:3*cLen/2);
                        aFilt = cFilt .* pFilt(cLen/2+1:3*cLen/2);
                    else,
                        cSig = sPyr{cLay}(c, :, id) .* pFilt(cLen+1:2*cLen) ./ iFilt2;
                        aFilt = cFilt .* pFilt(cLen+1:2*cLen) ./ iFilt2;
                    end
                    avg = sum(cSig) / sum(aFilt);
                    shifted = cSig - avg * aFilt;
                    antialias = exp( -(( 2*pi /cLen )^2) * ((movPyr{cLay}(r, c)^2)*fac*rho^2 + (nMov^2)*fac2*amp(id)^2*rhor^2) / 2);
                    %if cLay == nLay,
                    %    mMov = movPyr{cLay}(r, c)^2;
                    %    if c > 1,
                    %        mMov = max(mMov, movPyr{cLay}(r, c-1)^2);
                    %    end
                    %    if c < 2^cLay-1,
                    %        mMov = max(mMov, movPyr{cLay}(r, c+1)^2);
                    %    end
                    %    antialias = antialias .* exp( -(( 2*pi /cLen )^2) * (mMov*(dif .^2 - fac)) / 2);
                    %end
                    frq = fft(shifted) .* exp(2*pi*1i*amp(id)*nMov*dif/cLen) .* antialias;
                    shifted = real(ifft(frq));
                    shifted = shifted + avg * aFilt;

                    if cLay ~= 1,
                        if p * 2 - 1 == c,
                            sPyr{cLay - 1}(p, 1:cLen, id) = ...
                                sPyr{cLay - 1}(p, 1:cLen, id) + shifted .* cFilt ./ iFilt1;
                        elseif p * 2 == c,
                            sPyr{cLay - 1}(p, cLen/2+1:3*cLen/2, id) = ...
                                sPyr{cLay - 1}(p, cLen/2+1:3*cLen/2, id) + shifted .* cFilt;
                        else,
                            sPyr{cLay - 1}(p, cLen+1:2*cLen, id) = ...
                                sPyr{cLay - 1}(p, cLen+1:2*cLen, id) + shifted .* cFilt ./ iFilt2;
                        end
                    else,
                        synth(:, id) = (shifted ./ cFilt)';
                    end
                end
            end
        end
    end
end

function synth = motionAmp(sig, cmpSig, cmpRef, st, trim, mov, amp, rho, mmfilters, filt)
    rst = clip(round(st + mov), 1, length(sig) - trim + 1);
    res = st + mov - rst;

    s1 = cmpSig(st:st+trim-1) .* filt;
    s1 = s1 - (sum(s1) / sum(filt)) * filt;
    s2 = cmpRef(rst:rst+trim-1) .* filt;
    s2 = s2 - (sum(s2) / sum(filt)) * filt;
    so = sig(st:st+trim-1) .* filt;
    avg = (sum(so) / sum(filt));
    so = so - avg * filt;

    p1 = buildPyr(s1, mmfilters);
    p2 = buildPyr(s2, mmfilters);
    po = buildPyr(so, mmfilters);
    ang = {};
    for lay = 1:length(p1),
        ang = {ang{:} angle(p1{lay}./(p2{lay} + 1e-30))};
        if lay > 1 & lay < length(p1),
            for id = 1:length(ang{lay-1}),
                if abs(ang{lay-1}(id)) >= pi / 2,
                    ang{lay}(id*2-1) = ang{lay-1}(id) * 2;
                    ang{lay}(id*2) = ang{lay-1}(id) * 2;
                end
            end
        elseif lay == length(p1),
            for id = 1:length(ang{lay-1}),
                if abs(ang{lay-1}(id)) >= pi / 2,
                    ang{lay}(id) = ang{lay-1}(id) * 2;
                end
            end
        end
    end

    synth = zeros(trim, length(amp));

    for aid = 1:length(amp),
        a = amp(aid);
        pr = {};
        ang{lay} = ang{lay} + res * 2 * pi * (2^(aid-1))/ length(filt);
        for lay = 1:length(po),
            pr = {pr{:} po{lay} .* exp(i*ang{lay}*a - (( ang{lay}*rho ).^2) / 2)};
        end
        synth(:, aid) = reconPyr(pr, mmfilters) + avg * filt;
    end
end

function pyr = buildPyr(sig, filters)
    imdft = fft(sig); %DFT of image
    nFilts = size(filters, 1);
    pyr = {};
    for k = 1:nFilts
        curFilter = filters(k, :);
        tempDFT = curFilter.*imdft; % Transform domain
        indices = getIDX(curFilter);
        tempDFT = tempDFT(indices);
        
        curResult = ifft(tempDFT);
        pyr = {pyr{:} curResult};
    end
end

function rec = reconPyr(pyr, filters)
    n = size(filters, 2);
    imdft = zeros(1, n); %DFT of image
    nFilts = size(filters, 1);
    for k = 1:nFilts
        curFilter = filters(k, :);
        indices = getIDX(curFilter);
        tmpDFT = zeros(1, n);
        tmpDFT(indices) = fft(pyr{k});
        tmpDFT = tmpDFT .* curFilter;
        tmpDFT(n:-1:n/2+2) = conj(tmpDFT(2:n/2));
        imdft = imdft + tmpDFT;
    end
    rec = real(ifft(imdft));
end

function idx = getIDX(curFilter)
    m = max(find(curFilter > 0));
    l = length(curFilter);
    idx = zeros(1, l);
    idx(1:m) = 1;
    idx(l/2+1) = 1;
    idx(l:-1:l-m+2) = 1;
    idx = find(idx > 0);
end
