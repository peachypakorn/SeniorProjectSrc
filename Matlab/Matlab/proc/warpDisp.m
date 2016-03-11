function [synth, mask] = warpDisp(img, dspl, mul, maxGap, fillHole)
    nR = size(img, 1);
    nC = size(img, 2);
    nCol = size(img, 3);

    dspl1 = dspl(1:nR, 1:nC-1);
    dspl2 = dspl(1:nR, 2:nC);
    renderList = [];

    origLoc = repmat(1:nC-1, [nR, 1]);
    interpR = repmat((1:nR)', [1, nC-1]);

    for iG = -ceil(maxGap*abs(mul)):ceil(maxGap*abs(mul)),
        diffLoc = ceil(mul * dspl1) + iG;
        loc2 = 1 + mul * dspl2;
        loc1 = mul * dspl1;
        mask = (diffLoc - loc1) .* (loc2 - diffLoc) >= -eps;
        mask = and(mask, abs(dspl1 - dspl2) <= maxGap + eps);

        newLoc = round(origLoc + diffLoc);
        mask = and(mask, and(newLoc >= 1, newLoc <= nC));

        consDiff = loc2 - loc1;
        interpD = (dspl1 .* (loc2 - diffLoc) + dspl2 .* (diffLoc - loc1)) ./ consDiff;
        interpP = (origLoc .* (loc2 - diffLoc) + (origLoc+1) .* (diffLoc - loc1)) ./ consDiff;

        % choose the later point if equal
        maskEq = abs(loc2 - loc1) < eps;
        interpD(maskEq) = dspl2(maskEq);
        interpP(maskEq) = origLoc(maskEq) + 1;

        pushList = [(newLoc(mask)-1) * nR + interpR(mask),...
                    interpD(mask), interpP(mask), interpR(mask)];
        renderList = [renderList; pushList];
    end

    % No interpolation
    origLoc = repmat(1:nC, [nR, 1]);
    interpR = repmat((1:nR)', [1, nC]);

    diffLoc = mul * dspl;
    mask = abs(diffLoc - round(diffLoc)) < eps;

    newLoc = round(origLoc + diffLoc);
    mask = and(mask, and(newLoc >= 1, newLoc <= nC));

    interpD = dspl;
    interpP = origLoc;
    pushList = [(newLoc(mask)-1) * nR + interpR(mask),...
                interpD(mask), interpP(mask), interpR(mask)];
    renderList = [renderList; pushList];

    % Depth check
    renderList = sortrows(renderList);
    renderList = renderList(size(renderList, 1):-1:1, :);
    [~, ia, ~] = unique(renderList(:, 1), 'rows', 'legacy');
    renderList = renderList(ia, :);

    imgEx = zeros([nR, nC+1, nCol]);
    imgEx(:, 1:nC, :) = img;
    synth = zeros([nR, nC, nCol]);
    origLoc = (floor(renderList(:, 3)+1e-8)-1) * nR + renderList(:, 4);
    frac = renderList(:, 3) - floor(renderList(:, 3));

    for iCol = 0:nCol-1,
        synth(renderList(:, 1) + iCol * nR * nC) = ...
            imgEx(origLoc + iCol * nR * (nC+1)) .* (1 - frac) + ...
            imgEx(origLoc + nR + iCol * nR * (nC+1)) .* frac;
    end

    mask = ones(nR, nC);
    mask(renderList(:, 1)) = 0;

    if fillHole,
        synthDepth = zeros([nR, nC]);
        synthDepth(renderList(:, 1)) = renderList(:, 2);

        lv = zeros(nR, 1, nCol);
        ld = -1000*ones(nR, 1);
        imgD = -1000*ones(nR, nC);
        for iC = 1:nC,
            lmask = mask(:, iC);
            umask = lmask;
            nUp = length(find(umask));
            synth(find(umask), iC, :) = lv(find(umask), 1, :);
            imgD(find(umask), iC) = ld(find(umask), 1);
            lv(find(1-lmask), 1, :) = synth(find(1-lmask), iC, :);
            ld(find(1-lmask)) = synthDepth(find(1-lmask), iC);
        end
        lv = zeros(nR, 1, nCol);
        ld = -1000*ones(nR, 1);
        for iC = nC:-1:1,
            lmask = mask(:, iC);
            umask = and(lmask, ld > imgD(:, iC));
            nUp = length(find(umask));
            synth(find(umask), iC, :) = lv(find(umask), 1, :);
            imgD(find(umask), iC) = ld(find(umask), 1);
            lv(find(1-lmask), 1, :) = synth(find(1-lmask), iC, :);
            ld(find(1-lmask)) = synthDepth(find(1-lmask), iC);
        end
    end
end
%[xx, yy] = meshgrid(1:nC,1:nR);
%a = ceil(renderList(:, 1) / nR - eps);
%b = renderList(:, 1) - (a-1) * nR;
%c = imgEx(origLoc + iCol * nR * (nC+1)) .* (1 - frac) + ...
%    imgEx(origLoc + nR + iCol * nR * (nC+1)) .* frac;
%F = TriScatteredInterp(a,b,c);

%synth(:, :, iCol+1) = F(xx,yy);