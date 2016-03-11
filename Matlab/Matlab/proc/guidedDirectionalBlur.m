function [res, conf] = guidedDirectionalBlur(im, guid, sig_c, direction)
    nR = size(im, 1);
    nC = size(im, 2);
    nP = size(im, 3);
    res = zeros([nR, nC, nP]);
    conf = zeros([nR, nC]);

    if direction == 1,
        rC = nC-1:-1:1;
        pC = nC;
    else,
        rC = 2:nC;
        pC = 1;
    end

    div = 2*ones(nR, 1);
    div(1) = 1.5;
    div(nR) = 1.5;

    res(:, pC, :) = im(:, pC, :);
    for iC = rC,
        for dR = -1:1,
            rR = max(1, 1-dR):min(nR, nR-dR);
            upC = exp(-sum((guid(rR, iC, :) - guid(rR+dR, pC, :)).^2, 3) / (sig_c^2));
            res(rR, iC, :) = res(rR, iC, :) + ...
                             (im(rR, iC, :) .* repmat(1-upC, [1, 1, nP]) + ...
                             res(rR+dR, pC, :) .* repmat(upC, [1, 1, nP])) * (1-0.5*abs(dR)); 
            conf(rR, iC) = conf(rR, iC, :) + (1-upC + conf(rR+dR, pC) .* upC) * (1-0.5*abs(dR)); 
        end
        res(:, iC, :) = res(:, iC, :) ./ repmat(div, [1, 1, nP]);
        conf(:, iC, :) = conf(:, iC, :) ./ div;
        pC = iC;
    end
end