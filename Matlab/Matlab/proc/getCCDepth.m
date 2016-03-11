function [depth1, depth2] = getCCDepth(img1, img2, maxDepth)
    nR = size(img1, 1);
    nC = size(img1, 2);
    nCol = size(img1, 3);
    imgEx1 = zeros(nR, nC+2*maxDepth, nCol);
    imgEx2 = zeros(nR, nC+2*maxDepth, nCol);
    imgEx1(:, maxDepth + (1:nC), :) = img1;
    imgEx2(:, maxDepth + (1:nC), :) = img2;
    imgEx1(:, 1:maxDepth, :) = repmat(img1(:, 1, :), [1, maxDepth, 1]);
    imgEx1(:, nC+maxDepth+(1:maxDepth), :) = repmat(img1(:, nC, :), [1, maxDepth, 1]);
    imgEx2(:, 1:maxDepth, :) = repmat(img2(:, 1, :), [1, maxDepth, 1]);
    imgEx2(:, nC+maxDepth+(1:maxDepth), :) = repmat(img2(:, nC, :), [1, maxDepth, 1]);

    h = fspecial('gaussian', [15, 15], 3);
    mean1 = imfilter(imgEx1, h, 'symmetric');
    mean2 = imfilter(imgEx2, h, 'symmetric');
    sig1 = sqrt(imfilter(imgEx1 .^ 2, h, 'symmetric') - mean1 .^ 2);
    sig2 = sqrt(imfilter(imgEx2 .^ 2, h, 'symmetric') - mean2 .^ 2);

    depth1 = zeros(nR, nC);
    bestCC1 = -1e10 * ones(nR, nC);
    depth2 = zeros(nR, nC);
    bestCC2 = -1e10 * ones(nR, nC);
    for d = -maxDepth:maxDepth,
        mul = img1 .* imgEx2(:, maxDepth + d + (1:nC), :);
        cc = sum(imfilter(mul, h, 'symmetric') - ...
             mean1(:, maxDepth + (1:nC), :) .* mean2(:, maxDepth + d + (1:nC), :), 3) ./ ...
             sum(sig1(:, maxDepth + (1:nC), :) .* sig2(:, maxDepth + d + (1:nC), :), 3);
        mask = cc > bestCC1;
        depth1(mask) = d;
        bestCC1(mask) = cc(mask);
        mul = img2 .* imgEx1(:, maxDepth - d + (1:nC), :);
        cc = sum(imfilter(mul, h, 'symmetric') - ...
             mean2(:, maxDepth + (1:nC), :) .* mean1(:, maxDepth - d + (1:nC), :), 3) ./ ...
             sum(sig2(:, maxDepth + (1:nC), :) .* sig1(:, maxDepth - d + (1:nC), :), 3);
        mask = cc > bestCC2;
        depth2(mask) = d;
        bestCC2(mask) = cc(mask);
    end
end