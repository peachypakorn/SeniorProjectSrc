function outFilters = get1DCosFilters(dims, nFilters)
% Return nFilters + 1 cosine filters with equal spacing

ctr = ceil((dims+0.5)/2);
rad = abs(((1:dims)-ctr)./(dims/2));
mask = (1:dims)>ctr;

for ii = 1:nFilters+1,
    rad_i = rad - double(ii - 1) / nFilters;
    rad_i = rad_i * pi * nFilters / 2;
    rad_i = clip(rad_i, -pi / 2, pi / 2);
    outFilters{ii} = cos(rad_i).*mask;
end
outFilters{nFilters+2} = zeros(size(outFilters{1}));
outFilters{nFilters+2}(ctr) = 1.0;

