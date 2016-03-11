function dL = lagrDepthProp(imL, imR, dvs)
    color_gamma = 0.5;
    grad_gamma = 0.5;
    alpha = 0.1;
    sig_C = 0.02;
    min_col_dif = 0.03;
    min_grd_dif = 0.006;
    neighbor = 5;

    nH = size(imL, 1);
    nW = size(imL, 2);
    nLevel = ceil(log2(nW)) - 3;

    dL = zeros(nH, nW);
    minV = 1e6 * ones(nH, nW);

    nEx = max(abs(dvs));
    imREx = zeros(nH, nW+2*nEx, size(imL, 3));
    imREx(:, nEx+(1:nW), :) = imR;
    imREx(:, 1:nEx, :) = repmat(imR(:,1,:), [1, nEx, 1]);
    imREx(:, nW+nEx+(1:nEx), :) = repmat(imR(:,nW,:), [1, nEx, 1]);

    cv = zeros(nH, nW, numel(dvs));
    id = 1;
    for iL = 1:nLevel,
        for id = 1:numel(dvs),
            dv = dvs(id);
            imLC = imresize(imL, 2^(iL-nLevel));
            imRC = imresize(imREx(:, nEx+dv+(1:nW), :), 2^(iL-nLevel));
            if id == 1,
                cv = zeros([size(imLC, 1), size(imLC, 2), numel(dvs)]);
            end
            gL = colorGrad(imLC);
            gR = colorGrad(imRC);

            dif = (max(0, sum((imLC - imRC).^2, 3) - min_col_dif^2) .^ (color_gamma/2)) * alpha + ...
                  (max(0, sum((gL - gR).^2, 3) - min_grd_dif^2) .^ (grad_gamma/2)) * (1-alpha);
            dif = imgaussfilt(dif, 1);
            cv(:, :, id) = dif;
        end
        minV = 1e6*ones([size(cv, 1), size(cv, 2)]);
        minD = zeros([size(cv, 1), size(cv, 2)]);
        for id = 1:numel(dvs),
            cmp = cv(:,:,id);
            upd = cmp < minV;
            minV(upd) = cmp(upd);
            minD(upd) = dvs(id);
        end
        secV = 1e6*ones([size(cv, 1), size(cv, 2)]);
        for id = 1:numel(dvs),
            cmp = cv(:,:,id);
            upd = logical((cmp < secV) .* (abs(minD - dvs(id)) > neighbor-0.5));
            secV(upd) = cmp(upd);
        end
        conf = secV - minV;
        if iL == 1,
            dL = minD;
        else,
            alp = min(1, max(0, (conf-0.02)/0.03));
            dL = imresize(dL, size(minD)) .* (1-alp) + minD .* alp;
        end
        dL = guidedfilter_color(imLC, dL, [3,3], 1e-3);
        imshow(imresize((dL-dvs(1))/(dvs(numel(dvs))-dvs(1)), [nH, nW], 'nearest'));
        pause;
    end
end


function res = colorGrad(im)
    res = cat(3, gradient(im(:,:,1)), gradient(im(:,:,2)), gradient(im(:,:,3)));
end