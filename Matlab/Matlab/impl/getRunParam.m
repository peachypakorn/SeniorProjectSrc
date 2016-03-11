function runParam = getRunParam(testID, runID)
    runParam = struct();
    if runID == 1 || strcmp(runID, 'plain'),
        runParam.runName = 'plain';
        runParam.runID = 1;
        runParam.rho = 0.0;
        runParam.rho2D = 0.0;
        runParam.ht = 7;
        runParam.nA = 4;
        runParam.dA = 1;
        runParam.eurLagrResize = 0.25;
        runParam.lagrResize = 1;
        runParam.nPushPull = 0;
    elseif runID == 2 || strcmp(runID, 'aa1'),
        runParam.runName = 'aa1';
        runParam.runID = 2;
        runParam.rho = 1.0;
        runParam.rho2D = 0.5;
        runParam.ht = 7;
        runParam.nA = 4;
        runParam.dA = 1;
        runParam.eurLagrResize = 0.25;
        runParam.lagrResize = 1;
        runParam.nPushPull = 0;
    elseif runID == 3 || strcmp(runID, 'arctan3'),
        runParam.runName = 'arctan3';
        runParam.runID = 3;
        runParam.rho = 0.0;
        runParam.rho2D = 0.0;
        runParam.ht = 7;
        runParam.nA = 4;
        runParam.dA = 1;
        runParam.eurLagrResize = 0.25;
        runParam.lagrResize = 1;
        runParam.nPushPull = 0;
        runParam.remap = 'arctan';
        runParam.remapDepth = 3;
    elseif runID == 4 || strcmp(runID, 'exp-stitch-v2'),
        runParam.runName = 'exp-stitch-v2';
        runParam.runID = 4;
        runParam.rho = 0.5;
        runParam.rho2D = 0.0;
        runParam.ht = 7;
        runParam.nA = 4;
        runParam.dA = 1;
        runParam.eurLagrResize = 0.25;
        runParam.lagrResize = 1;
        runParam.nPushPull = 0;
        runParam.forcematch = true;
        runParam.sigWrap = 64;
        runParam.appWrap = 4;
    elseif runID == 5 || strcmp(runID, 'exp-stitch-noaa'),
        runParam.runName = 'exp-stitch-noaa';
        runParam.runID = 5;
        runParam.rho = 0.0;
        runParam.rho2D = 0.0;
        runParam.ht = 7;
        runParam.nA = 4;
        runParam.dA = 1;
        runParam.eurLagrResize = 0.25;
        runParam.lagrResize = 1;
        runParam.nPushPull = 0;
        runParam.forcematch = true;
        runParam.sigWrap = 64;
        runParam.appWrap = 4;
    elseif runID == 6 || strcmp(runID, 'stitch-lf'),
        runParam.runName = 'lf-r128-stitch-shear-arctan1';
        runParam.runID = 6;
        runParam.rho = 0.5;
        runParam.rho2D = 0.0;
        runParam.ht = 7;
        runParam.nA = 32;
        runParam.dA = 0.125;
        runParam.eurLagrResize = 0.25;
        runParam.lagrResize = 1;
        runParam.nPushPull = 0;
        runParam.forcematch = true;
        runParam.sigWrap = 64;
        runParam.appWrap = 32;
        runParam.recR = 128 + (0:1);

        runParam.expShear = true;
        runParam.shearLev = 4;
        runParam.shearMinLev = 1;
        runParam.shearGrid = [4, 4];

        runParam.remap = 'arctan';
        runParam.remapDepth = 1;
    elseif runID == 7 || strcmp(runID, 'full'),
        runParam.runName = 'full';
        runParam.runID = 7;
        runParam.rho = 0.5;
        runParam.rho2D = 0.25;
        runParam.ht = 9;
        runParam.nA = 4;
        runParam.dA = 1;
        runParam.eurLagrResize = 0.25;
        runParam.lagrResize = 1;
        runParam.nPushPull = 0;
        runParam.remap = 'arctan';
        runParam.remapDepth = 4;
        runParam.forcematch = true;
        runParam.sigWrap = 64;
        runParam.appWrap = 4;
    elseif runID == 8 || strcmp(runID, 'exp-sinc'),
        runParam.runName = 'exp-sinc';
        runParam.runID = 8;
        runParam.rho = 0.5;
        runParam.rho2D = 0.25;
        runParam.ht = 7;
        runParam.nA = 4;
        runParam.dA = 1;
        runParam.eurLagrResize = 0.25;
        runParam.lagrResize = 1;
        runParam.nPushPull = 0;
        runParam.expSinc = true;
    elseif runID == 9 || strcmp(runID, 'exp-shear'),
        runParam.runName = 'exp-shear';
        runParam.runID = 9;
        runParam.rho = 0.5;
        runParam.rho2D = 0.25;
        runParam.ht = 7;
        runParam.nA = 4;
        runParam.dA = 1;
        runParam.eurLagrResize = 0.25;
        runParam.lagrResize = 1;
        runParam.nPushPull = 0;
        runParam.expShear = true;
        runParam.shearLev = 4;
        runParam.shearMinLev = 1;
        runParam.shearGrid = [4, 4];
        runParam.forcematch = true;
        runParam.sigWrap = 64;
        runParam.appWrap = 4;
    elseif runID == 10 || strcmp(runID, 'full-shear'),
        runParam.runName = 'full-exp-shear';
        runParam.runID = 10;
        runParam.rho = 0.5;
        runParam.rho2D = 0.25;
        runParam.ht = 9;
        runParam.nA = 4;
        runParam.dA = 1;
        runParam.eurLagrResize = 0.25;
        runParam.lagrResize = 1;
        runParam.nPushPull = 0;
        runParam.remap = 'arctan';
        runParam.remapDepth = 4;
        runParam.expShear = true;
        runParam.shearLev = 6;
        runParam.shearMinLev = 2;
        runParam.shearGrid = [4, 4];
        runParam.forcematch = true;
        runParam.sigWrap = 64;
        runParam.appWrap = 4;
    else if runID == 11 || strcmp(runID, 'plain-v2'),
        runParam.runName = 'plain-v2';
        runParam.runID = 11;
        runParam.rho = 0;
        runParam.rho2D = 0.0;
        runParam.ht = 7;
        runParam.nA = 4;
        runParam.dA = 1;
        runParam.eurLagrResize = 0.25;
        runParam.lagrResize = 1;
        runParam.nPushPull = 0;
        runParam.expPrealign = true;
    elseif runID == 12 || strcmp(runID, 'full-v2'),
        runParam.runName = 'full-v2';
        runParam.runID = 12;
        runParam.rho = 0.5;
        runParam.rho2D = 0.25;
        runParam.ht = 8;
        runParam.nA = 4;
        runParam.dA = 1;
        runParam.eurLagrResize = 0.25;
        runParam.lagrResize = 1;
        runParam.nPushPull = 0;
        runParam.remap = 'arctan';
        runParam.remapDepth = 4;
    elseif runID == 13 || strcmp(runID, 'aa1-v2'),
        runParam.runName = 'aa1-v2';
        runParam.runID = 1;
        runParam.rho = 0.5;
        runParam.rho2D = 0.25;
        runParam.ht = 7;
        runParam.nA = 4;
        runParam.dA = 1;
        runParam.eurLagrResize = 0.25;
        runParam.lagrResize = 1;
        runParam.nPushPull = 0;
        runParam.expPrealign = true;
    elseif runID == 14 || strcmp(runID, 'halfhd-shear'),
        runParam.runName = 'halfhd-shear';
        runParam.runID = 14;
        runParam.rho = 0.5;
        runParam.rho2D = 0.25;
        runParam.ht = 8;
        runParam.nA = 4;
        runParam.dA = 1;
        runParam.eurLagrResize = 0.25;
        runParam.lagrResize = 1;
        runParam.nPushPull = 0;
        runParam.remap = 'arctan';
        runParam.remapDepth = 2;
        runParam.expShear = true;
        runParam.shearLev = 5;
        runParam.shearMinLev = 2;
        runParam.shearGrid = [4, 4];
        runParam.forcematch = true;
        runParam.sigWrap = 64;
        runParam.appWrap = 4;
    else,
        error('Run case not defined');
    end
end