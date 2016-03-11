function [dL, conf, minV] = bilateralDepth1S(imL, imR, dvs)
    border_diff = 0.3;
    color_gamma = 0.5;
    grad_gamma = 0.5;
    alpha = 0.1;
    border_weight = 1e-30;
    sig_C = 0.02;
    min_col_dif = 0.03;
    min_grd_dif = 0.003;

    nH = size(imL, 1);
    nW = size(imL, 2);
    nLevelW = ceil(log2(nW)) - 3;
    nLevelH = ceil(log2(nH)) - 3;

    dL = zeros(nH, nW);
    minV = 1e6 * ones(nH, nW);

    imL_g = rgb2gray(imL);
    imR_g = rgb2gray(imR);
    gdL = gradient(imL_g);
    gdR = gradient(imR_g);
    %close all;
    %figure;

    cv = zeros(nH, nW, numel(dvs));
    id = 1;
    for dv = dvs
        disp(dv);
        mask = zeros(nH, nW);
        yR = max(1, 1-dv):min(nW, nW-dv);
        mask(:, yR) = 1;
        dif = ones(nH, nW) * border_diff;
        dif(:, yR) = (max(0, sum((imL(:, yR, :) - imR(:, yR+dv, :)).^2, 3) - min_col_dif^2) .^ (color_gamma/2)) * alpha + ...
                      max(0, abs(gdL(:, yR) - gdR(:, yR+dv)) - min_grd_dif) .^ (grad_gamma) * (1-alpha);

        %dif(:, yR) = guidedfilter_color(imR(:, yR+dv, :), dif(:, yR), 9, 1e-4);
        %dif = multiscaleBilateral(dif, imL, (border_weight + mask), nLevelW, sig_C, 2, 1);
        %dif = multiscaleBilateral(dif, imL, (border_weight + mask), nLevelH, sig_C*2, 1, 2);
        %dif = multiscaleBilateral(dif, imL, (border_weight + mask), nLevelW, sig_C, 2, 1);
        dif = imgaussfilt(dif, 1);
        %dif = guidedfilter_color(imL, dif, [3,15], 1e-3);
        %dif = boxfilter(dif, 5);
        cv(:, :, id) = dif;
        id = id+1;

        if mod(dv, 10) == 0,
            %figure;
            imshow(dif);
            pause;
        end
    end

    for id = 1:numel(dvs),
        dv = dvs(id);
        if 0
            if id <= 3 || id > numel(dvs)-3
                continue;
            end
            dif = -min((cv(:,:,id+1:id+3) - repmat(cv(:,:,id), [1, 1, 3])) ./ repmat(reshape(1:3, [1,1,3]), [nH, nW, 1]), [], 3) +...
                   max((cv(:,:,id-3:id-1) - repmat(cv(:,:,id), [1, 1, 3])) ./ repmat(reshape(-3:-1, [1,1,3]), [nH, nW, 1]), [], 3);
            if mod(dv, 10) == 0,
                %figure;
                imshow(-dif*10+0.1);
                pause;
            end
        else
            dif = cv(:,:,id);
        end
        update = dif < minV;
        minV(update) = dif(update);
        dL(update) = dv;
    end

       % dL = multiscaleBilateral(dL, imL, (border_weight + mask), nLevelW, sig_C, 2, 1);
       % dL = multiscaleBilateral(dL, imL, (border_weight + mask), nLevelH, sig_C*2, 1, 2);
       % dL = multiscaleBilateral(dL, imL, (border_weight + mask), nLevelW, sig_C, 2, 1);
       %dL = guidedfilter_color(imL, dL, 11, 1e-4);

    secV = 1e6 * ones(nH, nW);
    for id = 1:numel(dvs),
        dif = cv(:, :, id);
        mask = min((dif < secV), (abs(dvs(id) - dL) > 4));
        secV(mask) = dif(mask);
    end
    conf = secV - minV;
end
