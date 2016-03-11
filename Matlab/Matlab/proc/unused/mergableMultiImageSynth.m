function [synth, depthP] = mergableMultiPyrSynth(img, imgref, amp, oc, rStrip, rStripDebug, trim, dscale, doffset)
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

    for lSt = 1:rStrip:size(img, 1),
        lEd = min(lSt + rStrip - 1, size(img, 1));
        movPyr = movePyrRec(cmpSig, cmpRef, oc, lSt, lEd, trim,...
            1, 1, size(img, 2), [1], 0);
        'pyrimid'  
        for r = lSt:lEd,
            for c = 1:3,
                [synth(r, :, c, :), dep] = mergableMultiPyrRec(img(r, :, c), movPyr, r-lSt+1, trim, ...
                        amp, 1, 1, 1, size(img, 2), ones([1 size(img, 2)]), 0);
                for l = 1:length(dep),
                    if l < length(dep),
                       depthP(r, :, :, length(dep)+1-l) = repmat(dep{l}' - dep{l+1}', 1, 3);
                    else,
                       depthP(r, :, :, length(dep)+1-l) = repmat(dep{l}', 1, 3);
                    end
                end
                depth(r, :, :) = repmat(dep{1}', 1, 3);
            end
            r
            if rStripDebug > 0,
                if mod(r, rStripDebug) == 0,
                    out_folder = sprintf('tmp_debug');
                    mkdir(out_folder);
                    for l = 1:length(dep),
                        imwrite(depthP(:, :, :, l) / min(dscale, size(img, 1) / 2^(l-1)) + 0.5, sprintf('%s/depth%03d.png', out_folder, l),'png');
                    end
                    imwrite(depth(:, :, :) / dscale + doffset, sprintf('%s/depth.png', out_folder),'png');
                    for m=1:length(amp),
                        %imwrite(outimg(:, :, :, m), sprintf('%s/%03d.png', out_folder, m),'png');
                        imwrite(repmat(synth(:, :, :, m), [1, 1, 1]), sprintf('%s/%03d.png', out_folder, m+2),'png');
                    end
                end
            end
        end
    end
end

function movPyr = movePyrRec(cmpSig, cmpRef, oc, lSt, lEd, trim, cLay, cN, cLen, cSt, cMov)
    oc = min(oc, cLen);
    sN = ((size(cmpSig, 2) / cLen) - 1) * oc + 1;
    spatialPyr = zeros([sN cLen]);
    spatialRef = zeros([sN cLen]);
    filt = 1:cLen;
    filt = cos((filt - (cLen+1)/2) * pi / cLen);
    %filt = filt - sum(filt) / cLen;
    %filt = ones(size(filt));

    sigma = cLen / 4;
    kernel = ceil(sigma * 3);
    if mod(kernel, 2) == 0,
        kernel = kernel + 1;
    end
    pSt = min(lSt - 1, (kernel - 1));
    pEd = min(size(cmpSig, 1) - lEd, (kernel - 1));
    nMov = zeros([kernel * 2 - 1 + lEd - lSt cN]);
    if numel(cMov) == 1,
        cMov = cMov * ones([kernel * 2 - 1 + lEd - lSt cN]);
    end
    for r = lSt-pSt:lEd+pEd,
        for idx = 1:sN,
            tmp = cmpSig(r, (1 + (idx-1) * cLen / oc):((idx-1+oc) * cLen / oc)) .* filt;
            off = sum(tmp) / sum(filt);
            tmp = tmp - off * filt;
            spatialPyr(idx, :) = fft(tmp);
            tmp = cmpRef(r, (1 + (idx-1) * cLen / oc):((idx-1+oc) * cLen / oc)) .* filt;
            off = sum(tmp) / sum(filt);
            tmp = tmp - off * filt;
            spatialRef(idx, :) = fft(tmp);
        end
        if cLen <= 4,
            samp = 1;
        elseif cLen <= 8,
            samp = 2;
        elseif cLen <= 16,
            samp = 4;
        else,
            samp = cLen/4;
        end
        p = 1:samp;
        for cI = 1:cN,
            st = cSt(cI);
            mov = cMov(r-lSt+1+kernel-1, cI);
            sId = round((st - 1) * 2 / cLen) * (oc / 2);
            rId = round((st + mov - 1) * oc / cLen);
            if rId < 0-oc/8 || rId >= sN+oc/8,
                nMov(r-lSt+1+kernel-1, cI) = 0;
            else,  
                rId = clip(rId, 0, sN-1);
                ang = angle(spatialPyr(sId+1, 1 + p) ./ (spatialRef(rId+1, 1 + p) + 1e-100));
                weight = abs(spatialPyr(sId+1, 1 + p) .* spatialRef(rId+1, 1 + p));
                for id = 4:samp,
                    angM = sum(ang(1:id-1).*p(1:id-1).*weight(1:id-1)) / sum(weight(1:id-1) .* (p(1:id-1).^2) + 1e-100);
                    ang(2:id) = ang(2:id) + round((angM * (2:id) - ang(2:id)) / (2*pi)) * 2 * pi;
                end
                angM = sum(ang .* p .* weight) / sum(weight .* (p.^2) + 1e-100);
                nMov(r-lSt+1+kernel-1, cI) = -mov + (rId - sId) * cLen / oc + angM * cLen / (2*pi);
            end
        end
    end
    for a = 1:kernel-1,
        if - a + lSt <= 0,
            nMov(kernel - a, :) = nMov(kernel + a - 1, :);
        end
        if lEd + a > size(cmpSig, 1),
            nMov(kernel + lEd - lSt + a, :) = nMov(kernel + lEd - lSt + 1 - a, :);
        end
    end
    nMov = imfilter(nMov + cMov, fspecial('gaussian', [floor(kernel/2), 1], sigma/2), 'symmetric') - cMov;
    %nMov = medfilt1(nMov + cMov, ceil(cLen / 4)) - cMov;
    
    if cLen <= trim,
        movPyr = {nMov(kernel:kernel+lEd-lSt, :)};
    else,
        rN = cN * 3;
        rSt = zeros(1, rN);
        rMov = zeros(kernel * 2 - 1 + lEd - lSt, rN);
        for cI = 1:cN,
            rSt(cI*3-2:cI*3) = cSt(cI) + (0:2) * (cLen/4);
            rMov(:, cI*3-2:cI*3) = repmat(cMov(:, cI) + nMov(:, cI), 1, 3);
        end
        nKernel = ceil((cLen / 8) * 3);
        if mod(nKernel, 2) == 0,
            nKernel = nKernel + 1;
        end
        rMov = rMov(kernel - nKernel+1:kernel + nKernel - 1 + lEd - lSt, :);
    cLen
        subPyr = movePyrRec(cmpSig, cmpRef, oc, lSt, lEd, trim, cLay+1, rN, cLen/2, rSt, rMov);
    cLen
        movPyr = {nMov(kernel:kernel+lEd-lSt, :), subPyr{:}};
    end
end

function [synth, dep] = mergableMultiPyrRec(sig, movPyr, r, trim, amp, cLay, cIdx, cSt, cLen, cF, mov)
    nA = length(amp);
    nMov = movPyr{cLay}(r, cIdx);
    if cLen <= trim,
        synth = repmat(sig(cSt:cSt+cLen-1)' .* cF', 1, length(amp));
        dep = {(mov + nMov) * cF};
    else,
        cF2 = 1:cLen/2;
        cF2 = cos((cF2 - (cLen/2+1)/2) * pi / (cLen/2));
        cF1 = zeros([1, cLen/2]);
        cF1(1:cLen/4) = 1;
        cF1(cLen/4+1:cLen/2) = cF2(cLen/4+1:cLen/2);
        cF3 = fliplr(cF1);
        [s1, d1] = mergableMultiPyrRec(sig, movPyr, r, trim, amp, ...
                cLay+1, cIdx * 3 - 2, cSt, cLen/2, cF(1:cLen/2).*cF1, mov+nMov);
        [s2, d2] = mergableMultiPyrRec(sig, movPyr, r, trim, amp, ...
                cLay+1, cIdx * 3 - 1, cSt+cLen/4, cLen/2, cF(1+cLen/4:3*cLen/4).*cF2, mov+nMov);
        [s3, d3] = mergableMultiPyrRec(sig, movPyr, r, trim, amp, ...
                cLay+1, cIdx * 3, cSt+cLen/2, cLen/2, cF(1+cLen/2:cLen).*cF3, mov+nMov);

        synth = zeros([cLen, nA]);
        synth(1:cLen/2, :) = synth(1:cLen/2, :) + s1 .* repmat(cF1', 1, nA);
        synth(cLen/4+1:3*cLen/4, :) = synth(cLen/4+1:3*cLen/4, :) + s2 .* repmat(cF2', 1, nA);
        synth(cLen/2+1:cLen, :) = synth(cLen/2+1:cLen, :) + s3 .* repmat(cF3', 1, nA);

        for lay = 1:length(d1),
            dep{lay} = zeros([1, cLen]);
            dep{lay}(1:cLen/2) = dep{lay}(1:cLen/2) + d1{lay} .* cF1;
            dep{lay}(cLen/4+1:3*cLen/4) = dep{lay}(cLen/4+1:3*cLen/4) + d2{lay} .* cF2;
            dep{lay}(cLen/2+1:cLen) = dep{lay}(cLen/2+1:cLen) + d3{lay} .* cF3;
        end
        dep{length(d1)+1} = cF * (nMov + mov);
    end

    dif = 0:cLen-1;
    dif(cLen/2+1) = 0;
    dif(cLen:-1:cLen/2+2) = -1:-1:-cLen/2+1;
    for id = 1:nA,
        avg = sum(synth(:, id)) / sum(cF);
        synth(:, id) = synth(:, id) - avg * cF';
        frq = fft(synth(:, id)) .* exp(2*pi*i*amp(id)*nMov*dif/cLen)';
        synth(:, id) = real(ifft(frq));
        synth(:, id) = synth(:, id) + avg * cF';
    end

        
end
        
