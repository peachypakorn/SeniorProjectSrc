function depth = eulaDepthProp(img1, img2)
    nR = size(img1, 1);
    nC = size(img1, 2);
    nLevel = ceil(log2(nC)) - 3;
    nExpand = 4;
    sig_c = 0.08;

    [filters, lowFilters] = getRadFilters(nC, nLevel);
    ratios = getDepthPhaseRatio(filters);

    pyr1 = cell(length(filters), 1);
    pyr2 = cell(length(filters), 1);
    low1 = cell(length(filters), 1);
    low2 = cell(length(filters), 1);
    pyrRef1 = cell(length(filters), 1);
    pyrRef2 = cell(length(filters), 1);
    lowRef1 = cell(length(filters), 1);
    lowRef2 = cell(length(filters), 1);
    imdft1 = fft(img1, [], 2);
    imdft2 = fft(img2, [], 2);
    for iL = 1:length(filters),
        pyr1{iL} = applyFilt(imdft1, filters{iL}, 1);
        pyr2{iL} = applyFilt(imdft2, filters{iL}, 1);
        low1{iL} = applyFilt(imdft1, lowFilters{iL}, 1);
        low2{iL} = applyFilt(imdft2, lowFilters{iL}, 1);
        pyrRef1{iL} = applyFilt(imdft1, filters{iL}, nExpand);
        pyrRef2{iL} = applyFilt(imdft2, filters{iL}, nExpand);
        lowRef1{iL} = applyFilt(imdft1, lowFilters{iL}, nExpand);
        lowRef2{iL} = applyFilt(imdft2, lowFilters{iL}, nExpand);
    end

    disp('pyramid construction done');
    depth = cell(length(filters), 1);
    depth{1} = zeros(nR, size(pyr1{1}, 2));
    for iL = 2:length(filters),
        disp(iL);
        imshow(imresize(colorRatio(pyr1{iL}), [nR, nC], 'nearest'));
        pause;
        depthUp = upsampleDepth(depth{iL-1}, colorRatio(pyr1{iL}), colorRatio(pyr1{iL-1}), 2, sig_c);
        nCL = size(pyr1{iL}, 2);
        [phsDif, weight, alig] = phaseDiff(pyr1{iL}, pyrRef2{iL}, ...
                low1{iL}, lowRef2{iL}, depthUp * nCL / nC, sig_c);
        imshow(imresize(abs(pyr1{iL}), [nR, nC], 'nearest')*10);
        pause;
        depth{iL} = alig * nC / nCL + phsDif * ratios{iL};
        imshow(imresize(depth{iL}, [nR, nC], 'nearest') / (2*pi* ratios{3}) + 0.5);
        pause;
        depth{iL} = guidedfilter_color(abs(pyr1{iL}), depth{iL}, [2^(length(filters) - iL)*2+1, 5], 1e-40);
        imshow(imresize(depth{iL}, [nR, nC], 'nearest') / (2*pi* ratios{3}) + 0.5);
        pause;
    end
end

function res = gradMod(im)
    res = mod(im - circshift(im, 1, 2), 2*pi);
end

function res = colorRatio(im)
    im = abs(im);
    res = im ./ sqrt(repmat(sum(im.^2, 3), [1,1,3]));
end