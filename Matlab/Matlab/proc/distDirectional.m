function dst = distDirectional(im, direction)
    nR = size(im, 1);
    nC = size(im, 2);
    nP = size(im, 3);
    dst = zeros(nR, nC);

    if direction == 1,
        rC = 1:nC;
    else,
        rC = nC:-1:1;
    end

    upd = zeros(nR, 1);
    pC = -1;
    for iC = rC,
        if pC ~= -1,
            upd1 = upd + sum(abs(im(:, iC, :) - im(:, pC, :)).^2, 3);
            if 0
                upd1(2:nR) = min(upd1(2:nR), upd(1:nR-1) + sum(abs(im(2:nR, iC, :) - im(1:nR-1, pC, :)), 3)*1.4);
                upd1(1:nR-1) = min(upd1(1:nR-1), upd(2:nR) + sum(abs(im(1:nR-1, iC, :) - im(2:nR, pC, :)), 3)*1.4);
            end
            upd = upd1;
        end
        upd1 = upd;
        if 1
            for iR = 2:nR,
                uv = upd(iR-1, 1) + sum(abs(im(iR-1, iC, :) - im(iR, iC, :)).^2);
                upd(iR, 1) = min(upd(iR, 1), uv);
            end
            for iR = nR-1:-1:1,
                uv = upd(iR+1, 1) + sum(abs(im(iR+1, iC, :) - im(iR, iC, :)).^2);
                upd(iR, 1) = min(upd(iR, 1), uv);
            end
        end
        %upd = max(upd, upd1-0.1);
        %upd = medfilt1(upd, 5);
        dst(:, iC) = upd;
        pC = iC;
    end
end