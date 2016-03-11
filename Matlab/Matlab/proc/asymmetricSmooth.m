function [filtered] = asymmetricSmooth(depth, fBack, fFront, sIni, sSmo, kFin, pow)
    kIni = floor(sIni * 1.5) * 2 + 1;
    kSmo = floor(sSmo * 1.5) * 2 + 1;
    hIni = fspecial('gaussian', [kIni, kIni], sIni);
    hSmo = fspecial('gaussian', [kSmo, kSmo], sSmo);

    depthBlur = imfilter(depth, hIni, 'symmetric');
    depthDiff = depth - depthBlur;
    depthDiff(depthDiff > 0) = depthDiff(depthDiff > 0) * fBack;
    depthDiff(depthDiff < 0) = -depthDiff(depthDiff < 0) * fFront;
    depthDiff = depthDiff + 1e-3;
    %imshow(depthDiff / 10);
    %figure;
    sigma = imfilter(depthDiff, hSmo, 'replicate');
    %imshow(sigma / 10);
    %figure;
    depth = exp(-depth * pow);
    filtered = adaptiveGaussian(depth, sqrt(sigma), kFin);
    filtered = log(filtered) * (-1/pow);
    %imshow(filtered / 40 + 0.5);
    %figure;
end