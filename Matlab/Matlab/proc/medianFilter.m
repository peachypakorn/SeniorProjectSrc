function [filtered] = medianFilter(image, weight, getShiftWeight, kernelR, kernelC, colProc)
    nR = size(image, 1);
    nC = size(image, 2);
    nCol = size(image, 3);

    hwR = floor(kernelR/2);
    hwC = floor(kernelC/2);
    nKE = (2*hwR+1)*(2*hwC+1);
    eIm = zeros(nR+2*hwR, nC+2*hwC, nCol);
    eWe = zeros(nR+2*hwR, nC+2*hwC);
    eIm(hwR+1:hwR+nR, hwC+1:hwC+nC, :) = image;
    eWe(hwR+1:hwR+nR, hwC+1:hwC+nC) = weight;

    filtered = zeros(nR, nC, nCol);
    for cSt = 1:colProc:nC,
        cEd = min(cSt+colProc-1, nC);
        cSize = cEd-cSt+1;
        idx = (1:cSize*nR)';
        for iCol = 1:nCol,
            collection = zeros([nKE*cSize*nR, 3]);
            iSt = 0;
            for dR = -hwR:hwR,
                for dC = -hwC:hwC,
                    val = eIm(hwR+dR+(1:nR), hwC+dC+(cSt:cEd), iCol);
                    we = eWe(hwR+dR+(1:nR), hwC+dC+(cSt:cEd)) .* ...
                         getShiftWeight(dR, dC, cSt, cEd);
                    val = reshape(val, [cSize*nR, 1]);
                    we = reshape(we, [cSize*nR, 1]);
                    newRows = [idx val we];
                    collection(iSt*cSize*nR + (1:cSize*nR), :) = newRows;
                    iSt = iSt+1;
                end
            end
            collection = sortrows(collection);
            collection = reshape(collection', [3, nKE, cSize*nR]);
            weSum = squeeze(cumsum(collection(3, :, :)));
            [~, rmin] = max(weSum > repmat(weSum(nKE, :), [nKE, 1]) * 0.5);
            val = collection(sub2ind(size(collection), 2*ones(1, cSize*nR), rmin, (1:cSize*nR)));
            filtered(:, cSt:cEd, iCol) = reshape(val, [nR, cSize]);
        end
    end
end