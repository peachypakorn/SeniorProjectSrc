function [outimg1, outimg2, nxtD] = prealign(img1, img2, maxD, preD, alpha)
    atc = zeros(maxD * 2 + 1);
    nC = size(img1, 2);
    for d = -maxD:maxD,
        atc(d+maxD+1) = sum(reshape((img1(:, maxD+1+d:nC-maxD+d, :) - img2(:, maxD+1:nC-maxD, :)).^2, [], 1));
    end
    [~, curD] = min(atc);
    curD = curD - maxD - 1;
    if abs(preD) <= maxD,
        nxtD = preD * alpha + curD * (1-alpha);
    else,
        nxtD = curD;
    end
    curD = round(nxtD);
    pL = floor((maxD - abs(curD)) / 2);
    pR = maxD - abs(curD) - pL;
    if curD <= 0,
        cL = 1 + pL;
        cR = nC - pR + curD;
    else,
        cL = 1 + pL + curD;
        cR = nC - pR;
    end
    outimg1 = img1(:, cL:cR, :);
    outimg2 = img2(:, cL-curD:cR-curD, :);
end