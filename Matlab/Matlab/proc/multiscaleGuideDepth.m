function [dL, conf, minV] = multiscaleGuideDepth(imL, imR, dvs)
    nH = size(imL, 1);
    nW = size(imL, 2);
    nLevel = ceil(log2(nW)) - 1;
    color_gamma = 0.25;

    cv = zeros(nH, nW, numel(dvs));
    dL = zeros(nH, nW);
    minV = 1e6 * ones(nH, nW);
    id = 1;
    for dv = dvs
        disp(dv);
        a = zeros(1,1,3);
        b = [0];
        for scale = 1:nLevel,
            imLS = imresize(imL, 2^(scale-nLevel));
            imRS = imresize(circshift(imR, dv, 2), 2^(scale-nLevel));
            aS = imresize(a, [size(imLS, 1), size(imLS, 2)]);
            bS = imresize(b, [size(imLS, 1), size(imLS, 2)]);
            dif = 2*(sum(aS.*imLS, 3) + bS) + (sum((imLS - imRS).^2, 3) .^ (color_gamma/2));

            [dif, a, b] = guidedfilter_color(imLS, dif, 2, 1e-4);
            if mod(dv, 10) == 0,
                %figure;
                imshow(dif*(2^(nLevel-scale))/1000);
                pause;
            end
        end
        cv(:, :, id) = dif;
        id = id+1;

        update = dif < minV;
        minV(update) = dif(update);
        dL(update) = dv;
    end

    secV = 1e6 * ones(nH, nW);
    for id = 1:numel(dvs),
        dif = cv(:, :, id);
        mask = min((dif < secV), (abs(dvs(id) - dL) > 4));
        secV(mask) = dif(mask);
    end
    conf = secV - minV;
end