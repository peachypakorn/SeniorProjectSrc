im1 = double(imread(sprintf('%sframe1.jpg', DATAPATH)))/255;
im2 = double(imread(sprintf('%sframe2.jpg', DATAPATH)))/255;

nR = 1080;
nC = 480;
im1 = imresize(im1, [nR, nC]);
im2 = imresize(im2, [nR, nC]);

kR = 7;
kC = 5;
maxDepth = 10;
[depth1, depth2, tmp1, tmp2, prevVals] = lagrDepthFast(im1, im2, -maxDepth, 1, 2*maxDepth+1, 2, 2, kR, kC, kR, kC, 1);
depth1 = maxDepth - depth1;
depth2 = maxDepth - depth2;

nEx = (512 - nC) / 2;
im1 = reflectEx(im1, nEx);
im2 = reflectEx(im2, nEx);
depth1 = reflectEx(depth1, nEx);
depth2 = reflectEx(depth2, nEx);

vId = 1;
testParam = struct();
testParam.maxDepth = maxDepth;
runParam = struct();
runParam.ht = 5;
runParam.nA = 8;
runParam.dA = 1;
runParam.rho = 0.3;
mem = 0;
phsCmpO = 0;
[synth, mem, phsCmpO, dump] = eulaOnLagrDumpResult(im1, im2, depth1, depth2, testParam, runParam, mem, vId, false, phsCmpO);
