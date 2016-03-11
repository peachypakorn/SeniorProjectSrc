function [synth] = mergableMultiMultiImageSynthFromDepth(img, depth, amp, oc, rStripDebug, trim)
    synth = zeros([ size(img) length(amp) ]);
    edepth = zeros([ size(img)]);

    cLen = size(img, 2);
    sdepth = {};
    while cLen >= trim,
        %sdepth = {sdepth{:}, imfilter(depth, fspecial('gaussian', [cLen, cLen], cLen/4), 'symmetric')};
        cLen
        div = 8;
        sigma = sqrt(((cLen/2)^2) / div);
        filtered = depth;
        for d = 1:4,
            filtered = max(filtered, imfilter(filtered, fspecial('gaussian', [1, cLen/2], sigma), 'symmetric'));
        end
        filtered = imfilter(filtered, fspecial('gaussian', [1, cLen], 2*sigma), 'symmetric');
        sdepth = {sdepth{:}, filtered};
        cLen = cLen / 2;
    end

    
    for r = 1:size(img, 1),
        for c = 1:3,
            [synth(r, :, c, :), edep] = mergableMultiPyrRec(img(r, :, c), sdepth, r, trim, ...
                    amp, 1, 1, 1, size(img, 2), ones([1 size(img, 2)]), 0);
            edepth(r, :, :) = repmat(edep', 1, 3);
        end
        r
        if rStripDebug > 0,
            if mod(r, rStripDebug) == 0,
                out_folder = sprintf('tmp_debug');
                mkdir(out_folder);
                imwrite((edepth / 255), sprintf('%s/depth.png', out_folder),'png');
                for m=1:length(amp),
                    %imwrite(outimg(:, :, :, m), sprintf('%s/%03d.png', out_folder, m),'png');
                    imwrite(repmat(synth(:, :, :, m), [1, 1, 1]), sprintf('%s/%03d.png', out_folder, m+2),'png');
                end
            end
        end
    end
end

function [synth, dep] = mergableMultiPyrRec(sig, sdepth, r, trim, amp, cLay, cIdx, cSt, cLen, cF, mov)
    nA = length(amp);
    
    nMov = max(sdepth{cLay}(r, cSt:cSt+cLen-1)) - mov;
    %nMov = (sdepth{cLay}(r, cSt+cLen/2-1) + sdepth{cLay}(r, cSt+cLen/2)) / 2 - mov;
    if cLen <= trim,
        synth = repmat(sig(cSt:cSt+cLen-1)' .* cF', 1, length(amp));
        dep = mov * cF;
    else,
        cF2 = 1:cLen/2;
        cF2 = cos((cF2 - (cLen/2+1)/2) * pi / (cLen/2));
        cF1 = zeros([1, cLen/2]);
        cF1(1:cLen/4) = 1;
        cF1(cLen/4+1:cLen/2) = cF2(cLen/4+1:cLen/2);
        cF3 = fliplr(cF1);
        [s1, d1] = mergableMultiPyrRec(sig, sdepth, r, trim, amp, ...
                cLay+1, cIdx * 3 - 2, cSt, cLen/2, cF(1:cLen/2).*cF1, mov+nMov);
        [s2, d2] = mergableMultiPyrRec(sig, sdepth, r, trim, amp, ...
                cLay+1, cIdx * 3 - 1, cSt+cLen/4, cLen/2, cF(1+cLen/4:3*cLen/4).*cF2, mov+nMov);
        [s3, d3] = mergableMultiPyrRec(sig, sdepth, r, trim, amp, ...
                cLay+1, cIdx * 3, cSt+cLen/2, cLen/2, cF(1+cLen/2:cLen).*cF3, mov+nMov);

        synth = zeros([cLen, nA]);
        synth(1:cLen/2, :) = synth(1:cLen/2, :) + s1 .* repmat(cF1', 1, nA);
        synth(cLen/4+1:3*cLen/4, :) = synth(cLen/4+1:3*cLen/4, :) + s2 .* repmat(cF2', 1, nA);
        synth(cLen/2+1:cLen, :) = synth(cLen/2+1:cLen, :) + s3 .* repmat(cF3', 1, nA);

        dep = zeros([1, cLen]);
        dep(1:cLen/2) = dep(1:cLen/2) + d1 .* cF1;
        dep(cLen/4+1:3*cLen/4) = dep(cLen/4+1:3*cLen/4) + d2 .* cF2;
        dep(cLen/2+1:cLen) = dep(cLen/2+1:cLen) + d3 .* cF3;
    end

    dif = 0:cLen-1;
    dif(cLen/2+1) = 0;
    dif(cLen:-1:cLen/2+2) = -1:-1:-cLen/2+1;
    for id = 1:nA,
        avg = sum(synth(:, id) ./ cF') / cLen;
        synth(:, id) = synth(:, id) - avg * cF';
        frq = fft(synth(:, id)) .* exp(2*pi*i*amp(id)*nMov*dif/cLen)';
        synth(:, id) = real(ifft(frq));
        synth(:, id) = synth(:, id) + avg * cF';
    end

        
end
        
