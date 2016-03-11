source_image_1 = double(imread('data/frame1.jpg')) / 255;
source_image_2 = double(imread('data/frame2.jpg')) / 255;

nLev = 5;
sigC = 0.08;
sigD = 5.5;
tau = 0.25;
lN = 11;
lN0 = 5;
lOmega = 10;

sizeExp = size(source_image_1);
sizeExp(1) = ceil(sizeExp(1) / (2^(nLev - 1))) * (2^(nLev - 1));
sizeExp(2) = ceil(sizeExp(2) / (2^(nLev - 1))) * (2^(nLev - 1));
exp_image_1 = zeros(sizeExp);
exp_image_2 = zeros(sizeExp);
exp_image_1(1:size(source_image_1, 1), 1:size(source_image_1, 2), :) = source_image_1;
exp_image_2(1:size(source_image_1, 1), 1:size(source_image_1, 2), :) = source_image_2;

[depthPyr] = simpleFlow1D(exp_image_1, exp_image_2, nLev, sigC, sigD, tau, lN, lOmega, lN0);

depthFactor = 32;
depthOffset = 0.5;
imwrite(depthPyr{1} / depthFactor + depthOffset, sprintf('%s/depth_l.png', 'result/tmp_debug'),'png');