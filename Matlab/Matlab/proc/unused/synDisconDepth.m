function img = synDisconDepth(h, v, dis, short, long)
% OUTPUT stereo image of hxv each view. left half has zero depth, and right 
%        half moved by dis. features white noise with period between short, long
% h, v, dis: int
    len = v + dis;
    if mod(len, 2) == 1,
        len = len + 1;
    end
    freq = zeros([1, len]);
    lo = clip(int32(len/long), 1, int32(len/2));
    hi = clip(int32(len/short), 1, int32(len/2));
    freq(lo+1:hi+1) = complex(rand([1, hi-lo+1])-0.5, rand([1, hi-lo+1])-0.5) * 0.25 * double(lo+hi+1.0) ./ double(lo+1:hi+1);
    freq(len-lo+1:-1:len-hi+1) = freq(len-lo+1:-1:len-hi+1) + conj(freq(lo+1:hi+1));
    signal = ifft(freq) * (double(len) / double(hi-lo+1)) + 0.5;

    img_1d = zeros([1, 2*v]);
    left = int32(v/2);
    img_1d(1:left) = signal(1:left);
    img_1d(left+1:v) = signal(left+1:v);
    %img_1d(v+1:v+left) = signal(1:left) + 0.03;
    %img_1d(v+left+1:2*v) = signal(left+1+dis:v+dis);
    img_1d(v+1:2*v) = signal(1+dis:v+dis);
    img = repmat(img_1d, [h, 1, 3]);
end
