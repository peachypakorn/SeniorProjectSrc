function [phsDif, weight, alig] = phaseDiff(wav, wavRef, low, lowRef, prealign, sig_c)
    nR = size(wav, 1);
    nC = size(wav, 2);
    nCRef = size(wavRef, 2);

    adim = size(wav);
    adim(1) = 1;
    adim(2) = 1;

    refIc = round((repmat(0:nC-1, [nR, 1]) + prealign) * (nCRef / nC));
    alig = refIc * (nC/nCRef) - repmat(0:nC-1, [nR, 1]);
    refIc = mod(refIc, nCRef) + 1;
    refIc = repmat(refIc, adim);
    [Ir, Ic, Icol] = ind2sub(size(wav), reshape(1:nR*nC*size(wav, 3), size(wav)));
    refInd = sub2ind(size(wavRef), Ir, refIc, Icol);
    phsDifCol = mod(angle(wav) - angle(wavRef(refInd))+pi, 2*pi)-pi;
    weightCol = min(abs(wav), abs(wavRef(refInd)));
    weight = sum(weightCol, 3);
    phsDif = sum(phsDifCol .* weightCol, 3) ./ weight;

    weightDif = exp(-sum((low - lowRef(refInd)).^2, 3) / (sig_c^2));
    weight = weight .* weightDif;
end