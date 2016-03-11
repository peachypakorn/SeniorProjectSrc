function [depthUp] = upsampleDepth(depth, lowUp, low, cMax, sig_c)
    nR = size(lowUp, 1);
    nC = size(lowUp, 2);
    nCRef = size(low, 2);

    depthUp = zeros(nR, nC);
    weight = zeros(nR, nC);

    for dx = -cMax:cMax,
        cId = mod(round((0:nC-1) * nCRef / nC) + dx, nCRef) + 1;
        wUpdate = exp(-sum((lowUp - low(:, cId, :)).^2, 3) / (sig_c^2));
        depthUp = depthUp + depth(:, cId, :) .* wUpdate;
        weight = weight + wUpdate;
    end

    depthUp = depthUp ./ weight;
end