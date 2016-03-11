function outimgs = lightfieldSubsample(ifiles, nA, rho)
    img = imread(ifiles{1});
    nR = size(img, 1);
    nC = size(img, 2);
    nCol = size(img, 3);
    nB = length(ifiles);
    outimgs = zeros([nR, nC, nCol, nA]);
    sig = nB * rho / nA + 1e-10;

    for iA = 1:nA,
        center = nB * (iA - 0.5) / nA + 0.5;
        dist = bsxfun(@min, abs((-nB+1:0) - center), ...
               bsxfun(@min, abs((1:nB) - center), abs((nB+1:2*nB) - center)));
        weight = exp(-dist .^ 2 / (2*sig^2));
        weight = weight / sum(weight);
        act_w = 0;
        for iB = 1:nB,
            if weight(iB) < 0.01,
                continue;
            else,
                img = double(imread(ifiles{iB})) / 255.0;
                act_w = act_w + weight(iB);
                outimgs(:, :, :, iA) = outimgs(:, :, :, iA) + img * weight(iB);
            end
        end
        outimgs(:, :, :, iA) = outimgs(:, :, :, iA) / act_w;
    end
end