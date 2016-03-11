function [synth] = mergableMultiPyrSynth(sig, sigref, amp, spatialPyr, spatialRef, spatialStart, oc)
    synth = mergableMultiPyrRec(sig, sigref, amp, spatialPyr, spatialRef, spatialStart, oc, ...
        1, 1, length(sig), ones([1, length(sig)]), 0);
end

function synth = mergableMultiPyrRec(sig, sigref, amp, spatialPyr, spatialRef, spatialStart, oc, ...
        cLay, cSt, cLen, cF, mov)
    sId = round((cSt - 1) * 2 / cLen) * (oc / 2);
    rId = round((cSt + mov - 1) * oc / cLen);
    if rId < 0 || rId >= spatialStart(3, cLay),
        nMov = 0;
    else,  
        if cLen <= 4,
            samp = 1;
        elseif cLen <= 8,
            samp = 2;
        elseif cLen <= 16,
            samp = 4;
        else,
            samp = 8;
        end
        p = 1:samp;
        ang = angle(spatialPyr(spatialStart(1, cLay) + spatialStart(2, cLay) * sId + 1 + p) ./ ...
              (spatialRef(spatialStart(1, cLay) + spatialStart(2, cLay) * rId + 1 + p) + 1e-100));
        for id = 4:samp,
            angM = sum(ang(1:id-1).*p(1:id-1)) / sum(p(1:id-1).^2);
            ang(id) = ang(id) + round((angM * id - ang(id)) / (2*pi)) * 2 * pi;
        end
        angM = sum(ang .* p) / sum(p.^2);
        %vs = spatialPyr(spatialStart(1, cLay) + spatialStart(2, cLay) * sId + 1+samp);
        %vr = spatialRef(spatialStart(1, cLay) + spatialStart(2, cLay) * rId + 1+samp);
        % angle(vr / (vs+1e-100))
        nMov = -mov + (rId - sId) * cLen / oc + angM * cLen / (2*pi);
    end

    nA = length(amp);
        
    if cLen == 2,
        synth = repmat(sig(cSt:cSt+1)' .* cF', 1, length(amp));
    else,
        cF2 = 1:cLen/2;
        cF2 = cos((cF2 - (cLen/2+1)/2) * pi / (cLen/2));
        cF1 = zeros([1, cLen/2]);
        cF1(1:cLen/4) = 1;
        cF1(cLen/4+1:cLen/2) = cF2(cLen/4+1:cLen/2);
        cF3 = fliplr(cF1);
        s1 = mergableMultiPyrRec(sig, sigref, amp, spatialPyr, spatialRef, spatialStart, oc, ...
                cLay+1, cSt, cLen/2, cF(1:cLen/2).*cF1, mov+nMov);
        s2 = mergableMultiPyrRec(sig, sigref, amp, spatialPyr, spatialRef, spatialStart, oc, ...
                cLay+1, cSt+cLen/4, cLen/2, cF(1+cLen/4:3*cLen/4).*cF2, mov+nMov);
        s3 = mergableMultiPyrRec(sig, sigref, amp, spatialPyr, spatialRef, spatialStart, oc, ...
                cLay+1, cSt+cLen/2, cLen/2, cF(1+cLen/2:cLen).*cF3, mov+nMov);

        synth = zeros([cLen, nA]);
        synth(1:cLen/2, :) = synth(1:cLen/2, :) + s1 .* repmat(cF1', 1, nA);
        synth(cLen/4+1:3*cLen/4, :) = synth(cLen/4+1:3*cLen/4, :) + s2 .* repmat(cF2', 1, nA);
        synth(cLen/2+1:cLen, :) = synth(cLen/2+1:cLen, :) + s3 .* repmat(cF3', 1, nA);
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
        
