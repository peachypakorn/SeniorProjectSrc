function [dL, dR, dL1, dL2, dR1, dR2, mL, mR] = lagrAdaptiveHorizontalBlur(imL, imR, dvs)
    border_diff = 0.3;
    color_gamma = 0.5;
    grad_gamma = 0.5;
    alpha = 0.1;
    min_col_dif = 0.00;
    min_grd_dif = 0.00;
    neighbor = 5;
    sig = 2;
    sig_c = 1e-30;%0.03;
    decay = 0.05;

    nH = size(imL, 1);
    nW = size(imL, 2);

    dL = zeros(nH, nW);
    dL1 = zeros(nH, nW);
    dL2 = zeros(nH, nW);
    minVL = 1e6 * ones(nH, nW);
    minVL1 = 1e6 * ones(nH, nW);
    minVL2 = 1e6 * ones(nH, nW);
    dR = zeros(nH, nW);
    dR1 = zeros(nH, nW);
    dR2 = zeros(nH, nW);
    minVR = 1e6 * ones(nH, nW);
    minVR1 = 1e6 * ones(nH, nW);
    minVR2 = 1e6 * ones(nH, nW);

    gL = colorGrad(imL);
    gR = colorGrad(imR);

    %cvL = zeros(nH, nW, numel(dvs));
    %cvR = zeros(nH, nW, numel(dvs));

    for id = 1:numel(dvs),
        dv = dvs(id);
        yR = max(1, 1-dv):min(nW, nW-dv);

        difL = ones(nH, nW) * border_diff;
        difL(:, yR) = (max(0, sum((imL(:, yR, :) - imR(:, yR+dv, :)).^2, 3) - min_col_dif^2) .^ (color_gamma/2)) * alpha + ...
                      (max(0, sum((gL(:, yR, :) - gR(:, yR+dv, :)).^2, 3) - min_grd_dif^2) .^ (grad_gamma/2)) * (1-alpha);
        difL = imgaussfilt(difL, sig);
        %difL = guidedfilter_color(imL, difL, 3, 1e-3);
        %[difL1, wL1] = horizontalBlur(difL, max(sqrt(sum(gL.^2, 3)) - min_grd_dif, 0), decay, 1);
        %[difL2, wL2] = horizontalBlur(difL, max(sqrt(sum(gL.^2, 3)) - min_grd_dif, 0), decay, -1);
        [difL1, wL1] = guidedDirectionalBlur(difL, imL, sig_c, 1);
        [difL2, wL2] = guidedDirectionalBlur(difL, imL, sig_c, -1);
        %cvL(:, :, id) = difL;

        update = difL1 < minVL1;
        minVL1(update) = difL1(update);
        dL1(update) = dv;
        update = difL2 < minVL2;
        minVL2(update) = difL2(update);
        dL2(update) = dv;

        difR = ones(nH, nW) * border_diff;
        difR(:, yR+dv) = (max(0, sum((imL(:, yR, :) - imR(:, yR+dv, :)).^2, 3) - min_col_dif^2) .^ (color_gamma/2)) * alpha + ...
                         (max(0, sum((gL(:, yR, :) - gR(:, yR+dv, :)).^2, 3) - min_grd_dif^2) .^ (grad_gamma/2)) * (1-alpha);
        difR = imgaussfilt(difR, sig);
        %difR = guidedfilter_color(imR, difR, 3, 1e-3);
        [difR1, wR1] = guidedDirectionalBlur(difR, imR, sig_c, 1);
        [difR2, wR2] = guidedDirectionalBlur(difR, imR, sig_c, -1);
        %cvR(:, :, id) = difR;

        update = difR1 < minVR1;
        minVR1(update) = difR1(update);
        dR1(update) = dv;
        update = difR2 < minVR2;
        minVR2(update) = difR2(update);
        dR2(update) = dv;
    end

    mask = dL1 < dL2;
    mask = min(mask, max(0, (wL1 - 0.5) / 0.3));
    mask = max(mask, min(1, 1 - (wL2 - 0.5) / 0.3));
    dL = round(dL1 .* mask + dL2 .* (1-mask));
    mask = dR1 < dR2;
    mask = min(mask, max(0, (wR1 - 0.5) / 0.3));
    mask = max(mask, min(1, 1 - (wR2 - 0.5) / 0.3));
    dR = round(dR1 .* mask + dR2 .* (1-mask));

    if 0
        secVL = 1e6 * ones(nH, nW);
        secVR = 1e6 * ones(nH, nW);
        for id = 1:numel(dvs),
            difL = cvL(:, :, id);
            mask = min((difL < secVL), (abs(dvs(id) - dL) > neighbor - 0.5));
            secVL(mask) = difL(mask);

            difR = cvL(:, :, id);
            mask = min((difR < secVR), (abs(dvs(id) - dR) > neighbor - 0.5));
            secVR(mask) = difR(mask);
        end
        confL = secVL - minVL;
        confR = secVR - minVR;
    end

    imshow((dL - min(dvs)) / (max(dvs) - min(dvs)));
    pause;
    imshow((dR - min(dvs)) / (max(dvs) - min(dvs)));
    pause;

    [mL, mR] = LRConsistencyCheck(dL, dR, 1);
    %confL = confL .* mL;
    %confR = confR .* mR;


    imshow(mL .* ((dL - min(dvs)) / (max(dvs) - min(dvs))));
    pause;
    imshow(mR .* ((dR - min(dvs)) / (max(dvs) - min(dvs))));
    pause;

end

function res = colorGrad(im)
    res = cat(3, gradient(im(:,:,1)), gradient(im(:,:,2)), gradient(im(:,:,3)));
end

function [res, weight] = horizontalBlur(val, conf, alpha, direction)
    nC = size(val, 2);
    nR = size(val, 1);
    nP = size(val, 3);
    if direction == 1,
        xC = 1:nC;
    else,
        xC = nC:-1:1;
    end

    res = zeros(size(val));
    weight = zeros(nR, nC);
    lV = zeros(nR, 1, nP);
    lW = zeros(nR, 1);
    for iC = xC,
        uW = exp(-conf(:, iC) / alpha);
        lV = lV .* repmat(uW, [1, 1, nP]) + val(:, iC, :) .* (1 - repmat(uW, [1, 1, nP]));
        lW = lW .* uW + 1 - uW;
        res(:, iC, :) = lV;
        weight(:, iC) = lW;

        %lV = imgaussfilt(lV, 1);
        %lW = imgaussfilt(lW, 1);
    end
    res = res ./ repmat(max(1e-30, weight), [1, 1, nP]);
end

function [mL, mR] = LRConsistencyCheck(dL, dR, th)
    [Y, X] = ind2sub(size(dL), reshape(1:numel(dL), size(dL)));

    XLR = min(size(dL, 2), max(1, X + dL));
    XRL = min(size(dR, 2), max(1, X - dR));

    mL = abs(dL - dR(sub2ind(size(dL), Y, XLR))) <= th;
    mR = abs(dR - dL(sub2ind(size(dL), Y, XRL))) <= th;
end