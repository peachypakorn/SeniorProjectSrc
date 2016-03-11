function [ans, iV, gV, wV] = multiscaleBilateral(im1, guid, weight, nLevel, sigC, dsW, dsH)
    gV = cell(nLevel, 1);
    wV = cell(nLevel, 1);
    iV = cell(nLevel, 1);
    gV{1} = guid;
    iV{1} = im1 .* weight;
    wV{1} = weight;
    for iL = 2:nLevel,
        [iV{iL}, gV{iL}, wV{iL}] = propagateUp(iV{iL-1}, gV{iL-1}, wV{iL-1}, sigC, dsW, dsH);
    end
    for iL = nLevel-1:-1:1,
        [iV{iL}, wV{iL}] = propagateDown(iV{iL+1}, gV{iL+1}, wV{iL+1}, ...
                                         iV{iL}, gV{iL}, wV{iL}, sigC);
    end
    ans = iV{1} ./ wV{1};
end

function [imU, guidU, wU] = propagateUp(im, guid, w, sigC, dsW, dsH)
    nH = size(im, 1);
    nW = size(im, 2);
    nHU = ceil(nH/dsH);
    nWU = ceil(nW/dsW);
    guidU = imresize(guid, [nHU, nWU]);
    wU = zeros(nHU,nWU);
    imU = zeros(nHU,nWU);
    for dx = -1:1,
        for dy = -1:1,
            guidSh = imresize(shiftArr(guidU, dx, dy), [nH, nW], 'nearest');
            wU = wU + shiftArr(imresize(w .* gaussValue(guidSh - guid, sigC), [nHU, nWU], 'box'), -dx, -dy);
            imU = imU + shiftArr(imresize(im .* gaussValue(guidSh - guid, sigC), [nHU, nWU], 'box'), -dx, -dy);
        end
    end
    wU = wU;
    imU = imU;
end

function [imD, wD] = propagateDown(imU, guidU, wU, imO, guidO, wO, sigC)
    nH = size(imO, 1);
    nW = size(imO, 2);
    nHU = size(imU, 1);
    nWU = size(imU, 2);
    wD = wO;
    imD = imO;
    for dx = -1:1,
        for dy = -1:1,
            guidSh = imresize(shiftArr(guidU, dx, dy), [nH, nW], 'nearest');
            wSh = imresize(shiftArr(wU, dx, dy), [nH, nW], 'nearest');
            imSh = imresize(shiftArr(imU, dx, dy), [nH, nW], 'nearest');
            wD = wD + wSh .* gaussValue(guidSh - guidO, sigC);
            imD = imD + imSh .* gaussValue(guidSh - guidO, sigC);
        end
    end
end

function imSh = shiftArr(im, dx, dy)
    imSh = zeros(size(im));
    xR = max(1, 1+dx):min(size(im, 1), size(im, 1)+dx);
    yR = max(1, 1+dy):min(size(im, 2), size(im, 2)+dy);
    imSh(xR, yR, :) = im(xR-dx, yR-dy, :);
end

function gv = gaussValue(dif, sigC)
    gv = exp(-sum(dif.^2, 3)/(sigC^2));
end