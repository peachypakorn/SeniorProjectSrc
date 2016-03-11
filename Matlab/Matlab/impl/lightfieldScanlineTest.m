%ipattern = '../bunnyOld1/%d.jpg';
ipattern = 'lf-result/bbb-bunny/lf-stitch-shear-v3/%d.jpg';
%ipattern = '~/Downloads/sup/data/bunny1/%d.jpg';
oprefix = 'lf-result/bbb-bunny/lf-stitch-shear-v3/sup';
%oprefix = '~/Downloads/sup/data/bunny1/sup';
mkdir(oprefix);

nStart = 0;
nI = 52;
ifiles = cell(nI, 1);
for iI = 1:nI,
    ifiles{iI} = sprintf(ipattern, iI+nStart-1);
end

nR = 166;
outimg = 0;
for iI = 1:nI,
    img = imread(ifiles{iI});
    if iI == 1,
        outimg = zeros(nI*2, size(img, 2), 3);
    end
    outimg(iI, :, :) = img(nR, :, :);
    outimg(iI+nI, :, :) = img(nR, :, :);
end

imwrite(outimg / 255, sprintf('%s/r%03d.jpg', oprefix, nR), 'jpg');