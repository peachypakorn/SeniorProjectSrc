function [outFilters, lowFilters] = getRadFilters(dims, nLev)
    [outFilters, lowFilters] = getRadFiltersImpl(dims, 2.^[-nLev+1:0], 2);
end

function [outFilters, lowFilters] = getRadFiltersImpl(dims, rVals, twidth)
    % Construct log_rad
    ctr = ceil((dims+0.5)/2);
    log_rad = abs(((1:dims)-ctr)./(dims/2));
    log_rad(ctr) = log_rad(ctr-1) / 4;

    N = max(size(rVals));

    [himaskPrev, lomask] = getRadialMaskPair(rVals(1),log_rad, twidth);
    mask = (1:dims)>ctr;
    outFilters = cell(N+1, 1);
    lowFilters = cell(N+1, 1);
    outFilters{1} = lomask;
    lowFilters{1} = outFilters{1};
    for k = 2:N
        [himask, lomask] = getRadialMaskPair(rVals(k),log_rad, twidth);
        outFilters{k} = lomask.*himaskPrev.*mask;
        lowFilters{k} = sqrt(lowFilters{k-1}.^2 + (lomask.*himaskPrev).^2);
        himaskPrev = himaskPrev.*himask;
    end
    outFilters{N+1} = himaskPrev.*mask;
    lowFilters{N+1} = sqrt(lowFilters{N}.^2 + himaskPrev.^2);
end

function [himask, lomask] = getRadialMaskPair( r, log_rad, twidth)
    log_rad  = log2(log_rad)-log2(r);
    twidth = min(twidth, -min(log_rad));

    himask = log_rad;
    himask = clip(himask, -twidth, 0 );
    himask = himask * pi/(2*twidth);
    himask = abs(cos(himask));
    lomask = sqrt(1-himask.^2);
end