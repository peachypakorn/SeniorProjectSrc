function waveletVis
    nC = 512;
    ht = 7;
    sizFilt = 2.^[0:-1:-ht];
    filters = get1DRadialFilters(nC, sizFilt, 1);
    nLev = max(size(filters));
    depPhsRatios = getDepthPhaseRatio(filters);
    pind = buildBandIndices(filters, 1);
    nF = pind(nLev, 1) + pind(nLev, 2) - 1;

    keyLev = 3;
    period = 12;
    sig1 = zeros(1, nC);
    sig2 = zeros(1, nC);
    for ii = 1:nC/2,
        if mod(ii-1, 2*period) < period,
            sig1(ii) = 1 + rand(1)/3;
        end
        if mod(ii-1-2, 2*period) < period,
            sig2(ii) = 1 + rand(1)/3;
        end
    end
    %sig1 = circshift(sig1, [0, 100]);
    %sig2 = circshift(sig2, [0, 100]);
    %sig2 = sig2 .* (1+(1:nC)/nC);

    pyr1 = buildSCFpyrGen1D(sig1, filters, 1);
    pyr2 = buildSCFpyrGen1D(sig2, filters, 1);

    ind = pyrBandIndices(pind, keyLev);
    phsCont = (nC / length(ind)) / (2^keyLev / (2*pi));
    maxS =  length(ind)/2;
    val = zeros(1, maxS);
    s1 = zeros(1, maxS);
    s2 = zeros(1, maxS);
    for iS = 1:maxS,
        nS = iS;
        tmpind = circshift(ind, [0, -floor(nS/2)]);
        s1(iS) = sum(pyr1(tmpind(1:nS)) .* exp(-i*(1:nS)'*phsCont) .* exp(-((1:nS) - (nS+1)/2)'.^2 / (2*(nS/2)^2)));
        s2(iS) = sum(pyr2(tmpind(1:nS)) .* exp(-i*(1:nS)'*phsCont) .* exp(-((1:nS) - (nS+1)/2)'.^2 / (2*(nS/2)^2)));
        val(iS) = angle(s2(iS)/s1(iS)) * (2^keyLev / (2*pi));
    end
    figure;
    hold on;
    plot(val);
    figure;
    hold on;
    plot(abs(s1), 'red');
    plot(abs(s2), 'black');

end

function outFilters = get1DRadialFilters(dims, rVals, twidth)
    % Construct log_rad
    ctr = ceil((dims+0.5)/2);
    log_rad = abs(((1:dims)-ctr)./(dims/2));
    log_rad(ctr) = log_rad(ctr-1) / 4;

    N = max(size(rVals));

    [himask, lomaskPrev] = getRadialMaskPair(rVals(1),log_rad, twidth);
    mask = (1:dims)>ctr;
    outFilters{1} = himask.*mask;
    for k = 2:N
        %figure;
        [himask, lomask] = getRadialMaskPair(rVals(k),log_rad, twidth);
        outFilters{k} = himask.*lomaskPrev.*mask;
        lomaskPrev = lomaskPrev.*lomask;
        %plot(real(ifft(ifftshift(outFilters{k}))))
    end
    outFilters{N+1} = lomaskPrev;
end

function [himask, lomask] = getRadialMaskPair( r, log_rad, twidth)
    if false, %XXX
        n = length(log_rad);
        ctr = ceil((n + 0.5)/2);
        sig = r * (n / 2);
        rang = max(1, ceil(ctr - 2*sig)):min(n, floor(ctr + 2*sig));
        lomask = zeros(size(log_rad));
        lomask(rang) = exp(- (rang - ctr).^2 / (2*sig^2));
        himask = sqrt(1 - lomask.^2);
        return;
    end
    log_rad  = log2(log_rad)-log2(r);

    himask = log_rad;
    himask = clip(himask, -twidth, 0 );
    himask = himask * pi/(2*twidth);
    himask = abs(cos(himask));
    lomask = sqrt(1-himask.^2);
end

function ratios = getDepthPhaseRatio(filters)
    nLev = max(size(filters));
    nC = length(filters{1});
    for k = 1:nLev,
        rotated = ifftshift(filters{k});
        ratios{k} = nC / (2 * pi * imag(sum(rotated .* (0:nC-1) .* 1i) / sum(rotated)));
    end
end

function pind = buildBandIndices(filters, expand)
    nLev = max(size(filters));
    pind = [];
    sumIdx = 1;
    for k = 1:nLev
        curFilter = filters{k};
        indices = getIDXFromFilter1D(curFilter, expand);
        pind = [pind; sumIdx, numel(indices)];
        sumIdx = sumIdx + numel(indices);
    end
end

function indices =  pyrBandIndices(pind,band)
    indices = pind(band, 1):pind(band, 1) + pind(band, 2) - 1;
end

function pyr = buildSCFpyrGen1D(im, filters, expand)
    % Return pyramid in the usual format of a stack of column vectors
    imdft = fftshift(fft(im)); %DFT of image
    nLev = max(size(filters));
    n = length(im);
    pyr = [];
    for k = 1:nLev
        curFilter = filters{k};
        tempDFT = curFilter.*imdft; % Transform domain

        indices = getIDXFromFilter1D(curFilter, expand);
        nUnex = length(getIDXFromFilter1D(curFilter, 1));
        % shift to avoid cross-boundary wavelet
        tempDFT = tempDFT .* exp(2*pi*1i*((0:n-1) - floor(n/2))/nUnex);
        tempDFT = tempDFT(indices);
        curResult = ifft(ifftshift(tempDFT)) * length(indices) / nUnex;
        curResult = circshift(curResult, [0, floor(length(indices)/(2*nUnex))]);
        pyr = [pyr; curResult(:)];
    end
end

function res = reconSCFpyrGen1D(pyr, pind, filters, minLev, maxLev)
    nLev = max(size(filters));
    n = length(filters{1});
    imdft = zeros(1, n); %DFT of image
    for k = minLev:maxLev
        curFilter = filters{k};
        tempDFT = fftshift(fft(2*real(pyr(pyrBandIndices(pind,k))'))); % Transform domain
        if k < nLev,
            curFilter = curFilter + [0, rot90(curFilter(2:n),2)];
        else
            tempDFT = tempDFT / 2;
        end
        indices = getIDXFromFilter1D(curFilter, 1);
        imdft(indices) = imdft(indices) + tempDFT.*curFilter(indices)...
                .*exp(-2*pi*1i*(indices-1-floor(n/2))/length(indices));
    end
    res =  real(ifft(ifftshift(imdft)));
end

function [ out, lenO ] = getIDXFromFilter1D(filt, expand)
    aboveZero = filt>1e-10;
    aboveZero = logical(aboveZero + fliplr(aboveZero));
    n = length(filt);

    idx1 = 1:length(filt);
    idx1 = idx1(aboveZero);
    lenO = max(idx1) - min(idx1) + 1;
    len = lenO;

    fac = 1;
    while fac * 2 <= expand && len * fac * 2 <= n,
        fac = fac * 2;
    end
    len = len * fac;
    ctr = ceil((n + 1) / 2);
    if mod(len, 2) == 0,
        out = ctr - len/2:ctr + len/2-1;
    else,
        out = ctr - (len-1)/2:ctr + (len-1)/2;
    end
end
