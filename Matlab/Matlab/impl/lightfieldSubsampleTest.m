%ipattern = '../bunnyOld1/%d.jpg';
ipattern = 'lf-result/train/lf-stitch-shear-cut/%d.jpg';
%ipattern = '~/Downloads/sup/data/train1/%d.jpg';
oprefix = 'result/lf8/train/stitch-shear-cut';
mkdir(oprefix);

nStart = 0;
nI = 27;
ifiles = cell(nI, 1);
for iI = 1:nI,
    ifiles{iI} = sprintf(ipattern, iI+nStart-1);
end

nA = 8;
rho = 0.5;
outimgs = lightfieldSubsample(ifiles, nA, rho);
for iA = 1:nA,
    imwrite(outimgs(:, :, :, iA), sprintf('%s/001_%03d.jpg', oprefix, iA), 'jpg');
end
