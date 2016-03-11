function [depthPyr] = simpleFlow1D(imgSig, imgRef, nLev, sigC, sigD, tau, lN, lOmega, lN0)
    depthPyr = cell(nLev, 1);
    irregular = 0;
    hN = (lN - 1) / 2;
    hN0 = (lN0 - 1) / 2;
    t0 = cputime();
    for iLev = nLev:-1:1,
        img1 = imresize(imgSig, 2^(-iLev+1));
        img2 = imresize(imgRef, 2^(-iLev+1));
        nR = size(img1, 1);
        nC = size(img1, 2);
        if iLev == nLev,
            depthPyr{iLev} = zeros(nR, nC);
            irregular = nLev * ones(nR, nC);
        else,
            depthPyr{iLev} = imresize(depthPyr{iLev+1}, 2) * 2;
            irregular = imresize(irregular, 2, 'nearest');
        end
        for iR = 1:nR,
            for iC = 1:nC,
                gr = 2^(irregular(iR, iC) - iLev);
                if mod(iR-1, gr) == 0 & mod(iC-1, gr) == 0,
                    depthPyr{iLev}(iR, iC) = findDepth(img1, img2, iR, iC, ...
                            depthPyr{iLev}(iR, iC), sigC, sigD, hN, lOmega);
                end
            end
        end
        for jLev = nLev:-1:iLev+1,
            gr = 2^(jLev - iLev);
            for iR = 1:gr:nR,
                for iC = 1:gr:nC,
                    if irregular(iR, iC) == jLev,
                        if iR + gr <= nR,
                            iR2 = iR + gr;
                        else,
                            iR2 = iR + gr - 1;
                            depthPyr{iLev}(iR2, iC) = findDepth(img1, img2, iR2, iC, ...
                                    depthPyr{iLev}(iR2, iC), sigC, sigD, hN, lOmega);
                        end
                        if iC + gr <= nC,
                            iC2 = iC + gr;
                        else,
                            iC2 = iC + gr - 1;
                            depthPyr{iLev}(iR, iC2) = findDepth(img1, img2, iR, iC2, ...
                                    depthPyr{iLev}(iR, iC2), sigC, sigD, hN, lOmega);
                            depthPyr{iLev}(iR2, iC2) = findDepth(img1, img2, iR2, iC2, ...
                                    depthPyr{iLev}(iR2, iC2), sigC, sigD, hN, lOmega);
                        end
                        od = 0:gr-1;
                        rdR = (iR2 - iR) - od;
                        rdC = (iC2 - iC) - od;
                        depthPyr{iLev}(iR:iR+gr-1, iC:iC+gr-1) = ( ...
                                depthPyr{iLev}(iR, iC)*rdR'*rdC + ...
                                depthPyr{iLev}(iR, iC2)*rdR'*od + ...
                                depthPyr{iLev}(iR2, iC)*od'*rdC + ...
                                depthPyr{iLev}(iR2, iC2)*od'*od) / ((iR2 - iR) * (iC2 - iC));

                    end
                end
            end
        end
        minMap = zeros(nR, nC);
        maxMap = zeros(nR, nC);
        for iR = 1:nR,
            for iC = 1:nC,
                minMap(iR, iC) = min(depthPyr{iLev}(iR, max(iC-hN0, 1):min(iC+hN0, nC)));
                maxMap(iR, iC) = max(depthPyr{iLev}(iR, max(iC-hN0, 1):min(iC+hN0, nC)));
            end
        end
        minMap2 = zeros(nR, nC);
        maxMap2 = zeros(nR, nC);
        for iR = 1:nR,
            for iC = 1:nC,
                minMap2(iR, iC) = min(minMap(max(iR-hN0, 1):min(iR+hN0, nR), iC));
                maxMap2(iR, iC) = max(maxMap(max(iR-hN0, 1):min(iR+hN0, nR), iC));
            end
        end
        for iR = 1:nR,
            for iC = 1:nC,
                if minMap2(iR, iC) < depthPyr{iLev}(iR, iC) - tau | ...
                        maxMap2(iR, iC) > depthPyr{iLev}(iR, iC) + tau,
                    while irregular(iR, iC) >= iLev,
                        gr = 2^(irregular(iR, iC)-iLev);
                        sR = floor((iR - 1) / gr) * gr;
                        sC = floor((iC - 1) / gr) * gr;
                        irregular(sR+1:sR+gr, sC+1:sC+gr) = irregular(sR+1:sR+gr, sC+1:sC+gr) - 1;
                    end
                end
            end
        end

        iLev
        cputime() - t0
        t0 = cputime();
        imwrite(irregular / nLev, sprintf('%s/irr-%d.png', 'result/tmp_debug', iLev),'png');
        imwrite(depthPyr{iLev}*(2^(iLev - 1)) / 32 + 0.5, sprintf('%s/depth-%d.png', 'result/tmp_debug', iLev),'png');
    end
end

function depth = findDepth(img1, img2, iR, iC, depth0, sigC, sigD, hN, lOmega)
    weight = zeros(2*hN+1, 2*hN+1);
    nR = size(img1, 1);
    nC = size(img1, 2);
    for iX = -hN:hN,
        for iY = -hN:hN,
            if iR + iX > 0 & iR + iX <= nR & iC + iY > 0 & iC + iY <= nC,
                weight(hN+1+iX, hN+1+iY) = exp(-(iX^2+iY^2)/(2*sigD^2) - ...
                        sum((img1(iR+iX, iC+iY, :) - img1(iR, iC, :)).^2)/(2*sigC^2));
            end
        end
    end
    arr = 1e3 * ones(1, lOmega);
    for iD = floor(depth0) - lOmega/2+1: floor(depth0) + lOmega/2,
        if iC + iD <= 0 | iC + iD > nC,
            continue;
        end
        sumW = 0;
        sumE = 0;
        for iY = -hN:hN,
            if iC + iD + iY > 0 & iC + iD + iY <= nC,
                for iX = -hN:hN,
                    if weight(hN+1+iX, hN+1+iY) == 0,
                        continue;
                    end
                    sumW = sumW + weight(hN+1+iX, hN+1+iY);
                    sumE = sumE + weight(hN+1+iX, hN+1+iY) * ...
                            sum((img2(iR+iX, iC+iY+iD, :) - img1(iR+iX, iC+iY, :)).^2);
                end
            end
        end
        if sumW > 0,
            arr(iD-(floor(depth0) - lOmega/2)) = sumE / sumW;
        end
    end
    [tmp, depth] = min(arr);
    %if (depth == 1 | depth == lOmega),
    %    [iR, iC]
    %    arr
    %end
    if depth + 1 <= lOmega & depth - 1 > 0,
        if arr(depth+1) < 1e9 & arr(depth - 1) < 1e9,
            y1 = arr(depth-1);
            y2 = arr(depth);
            y3 = arr(depth+1);
            a = ((y3 - y2) + (y1 - y2)) / 2;
            depth = (-(y3-y2) + a*(2*depth+1)) / (2*a);
        end
    end
    depth = depth + (floor(depth0) - lOmega/2);
end