%source_image_1 = double(imread('data/frame1.jpg')) / 255;
%source_image_2 = double(imread('data/frame2.jpg')) / 255;
%source_image_1 = double(imread('data/bbb_69_left_062.jpg')) / 255;
%source_image_2 = double(imread('data/bbb_69_right_062.jpg')) / 255;
source_image_1 = double(imread('data/bbb_252p7_right.jpg')) / 255;
source_image_2 = double(imread('data/bbb_252p7_left.jpg')) / 255;
%source_image_1 = double(imread('data/bbb_305_right_228.jpg')) / 255;
%source_image_2 = double(imread('data/bbb_305_left_228.jpg')) / 255;
%source_image = double(imread('data/video_2_undist_1.jpg')) / 255;
%source_image = double(imread('data/transparency_undist.jpg')) / 255;
%source_image_1 = source_image(:, size(source_image, 2)/2+1:size(source_image, 2), :);
%source_image_2 = source_image(:, 1:size(source_image, 2)/2, :);

%source_image_1 = double(imread('data/transparency_right.jpg')) / 255;
%source_image_2 = double(imread('data/transparency_left.jpg')) / 255;
%source_image_1 = double(imread('data/check9_1.png')) / 255;
%source_image_2 = double(imread('data/check9_2.png')) / 255;
%source_image_1 = double(imread('data/art_1.png')) / 255;
%source_image_2 = double(imread('data/art_2.png')) / 255;
%source_image = double(imread('data/transparency_undist.jpg')) / 255;
%source_image_1 = source_image(:, 1:size(source_image, 2)/2, :);
%source_image_2 = source_image(:, size(source_image, 2)/2+1:size(source_image, 2), :);
source_image_1 = imresize(source_image_1, 1/3);
source_image_2 = imresize(source_image_2, 1/3);
%source_image_1 = imresize(source_image_1, 1/2);
%source_image_2 = imresize(source_image_2, 1/2);


source_image_1 = source_image_1(:, :, :);
source_image_2 = source_image_2(:, :, :);
%depthL = 0;
%depthR = 0;
depthL = round(double(imread('data/bbb-252p7_depthL.jpg'))-128);
depthR = round(double(imread('data/bbb-252p7_depthR.jpg'))-128);
%depthL = round(double(imread('data/art_d.png')));
%depthL = fillPixelsReference(source_image_1, depthL, 0.1, 9, 19, -1);
%depthL = -depthL/4;
%depthR = warpImage(depthL, depthL, depthL);

rho = 0.0;
ht = 6;
nA = 4;
dA = 1;

outimg = zeros([ size(source_image_1) 2*nA ]);
t0 = cputime();
%[outimg, depthL, depthR] = lagrangian(source_image_1, source_image_2, nA, dA, 1, 10, 20, depthL, depthR);
%[outimg, depthL, depthR] = waveLagr(source_image_1, source_image_2, rho, ht, nA, dA, 1, 10, 20, 0, 0, 0, depthL, depthR);
[outimg, depthL, depthR] = laplLagr(source_image_1, source_image_2, ht, nA, dA, 1, 20, 40, depthL, depthR);
%outimg = eulaOnLagr(source_image_1, source_image_2, rho, ht, nA, dA, 1, 20, 40, 0, 0, 0);
cputime() - t0

out_folder = sprintf('result/lagrangian/bbb-252p7/laplwarp');
%out_folder = sprintf('result/tmp/test-phase');
out_folder = sprintf('%s_%.2f', out_folder, rho);
mkdir(out_folder);
for m=1:2*nA,
    imwrite(outimg(:, :, :, m), sprintf('%s/%03d.jpg', out_folder, m),'jpg');
end
imwrite((depthL + 128)/255, sprintf('%s/depthL.jpg', out_folder),'jpg');
imwrite((depthR + 128)/255, sprintf('%s/depthR.jpg', out_folder),'jpg');