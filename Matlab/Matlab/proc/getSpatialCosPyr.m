function [pyr, start] = getSpatialCosPyr(sig, oc)
% only support power of two, oc is also power of two (at least 2 for overlap)
    % calculate array sizes
    n = length(sig);
    nl = log2(n);
    np = 0;
    for l = 1:nl,
        fl = n / (2^(l-1));
        ns = (2^(l-1) - 1) * oc + 1;
        np = np + ns * fl;
    end
    pyr = zeros([1 np]);
    start = zeros([3 nl+1]);
    % actual generation
    ed = 0;
    for l = 1:nl,
        fl = n / (2^(l-1));
        ns = (2^(l-1) - 1) * oc + 1;
        start(1, l) = ed;
        start(2, l) = fl;
        start(3, l) = ns;
        filt = 1:fl;
        filt = cos((filt - (fl+1)/2) * pi / fl);
        filt = filt - sum(filt) / fl;
        filt = ones(size(filt));
        for idx = 1:ns,
            tmp = sig((1 + (idx-1) * fl / oc):((idx-1+oc) * fl / oc)) .* filt;
            pyr(ed+1:ed+fl) = fft(tmp);
            ed = ed + fl;
        end
    end
end
            
        
    
