%ipattern = '../bunnyOld1/%d.jpg';
ipattern = '~/Downloads/sup/data/bunny1/%d.jpg';
%oprefix = 'lf-result/bbb-bunny/lf-shear-stitch-v3';
oprefix = 'lf-result/bbb-bunny/lf-stitch-shear-v4';
mkdir(oprefix);

nI = 100;
ifiles = cell(nI, 1);
for iI = 1:nI,
    ifiles{iI} = sprintf(ipattern, iI-1);
end

runParam = struct();
runParam.ht = 7;

runParam.sigWrap = 64;
%runParam.sigWrap = 64 * 8;
runParam.appWrap = 24;
runParam.reduWrap = 24;

runParam.expShear = true;
runParam.shearLev = 4;
runParam.shearMinLev = 2;
runParam.shearGrid = [4, 4];
%runParam.maxShear = 100;

t0 = cputime();
phsCmpO = lightfieldStitchV2(ifiles, oprefix, runParam);
cputime() - t0;
