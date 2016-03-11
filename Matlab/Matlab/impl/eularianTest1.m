
source_image_1 = double(imread('data/frame1.jpg')) / 255;
source_image_2 = double(imread('data/frame2.jpg')) / 255;
%source_image_1 = double(imread('data/bbb_69_left_062.jpg')) / 255;
%source_image_2 = double(imread('data/bbb_69_right_062.jpg')) / 255;
%source_image_1 = double(imread('data/bbb_252p7_right.jpg')) / 255;
%source_image_2 = double(imread('data/bbb_252p7_left.jpg')) / 255;
%source_image_1 = double(imread('data/bbb_305_right_228.jpg')) / 255;
%source_image_2 = double(imread('data/bbb_305_left_228.jpg')) / 255;
%source_image = double(imread('data/video_2_undist_1.jpg')) / 255;
%source_image = double(imread('data/transparency_undist.jpg')) / 255;
%source_image_1 = source_image(:, size(source_image, 2)/2+1:size(source_image, 2), :);
%source_image_2 = source_image(:, 1:size(source_image, 2)/2, :);

%source_image_1 = double(imread('data/transparency_right.jpg')) / 255;
%source_image_2 = double(imread('data/transparency_left.jpg')) / 255;
%source_image_1 = double(imread('data/window1_1.jpg')) / 255;
%source_image_2 = double(imread('data/window1_2.jpg')) / 255;
%source_image_1 = double(imread('data/window1_1.jpg')) / 255;
%source_image_2 = double(imread('data/window1_2.jpg')) / 255;
%source_image = double(imread('data/bicycle.jpg')) / 255;
%source_image_1 = source_image(:, 1:size(source_image, 2)/2, :);
%source_image_2 = source_image(:, size(source_image, 2)/2+1:size(source_image, 2), :);
%source_image_1 = imresize(source_image_1, 2/3);
%source_image_2 = imresize(source_image_2, 2/3);

%source_image_1 = imresize(source_image_1, 1/2);
%source_image_2 = imresize(source_image_2, 1/2);

rho = 1.0;
ht = 6;
nA = 4;
dA = 3;

outimg = zeros([ size(source_image_1) 2*nA ]);
t0 = cputime();
%outimg = eularian1D(source_image_1, source_image_2, rho, ht, nA, dA, 0, 0);
outimg = eulaOnLagr(source_image_1, source_image_2, rho, ht, nA, dA, 0.25, 10, 20, 3, 0);
cputime() - t0

out_folder = sprintf('result/eularian-1D/bbb-basic/forcematch-wrap-amp3-arctan3');
%out_folder = sprintf('result/tmp/test-phase');
out_folder = sprintf('%s_%.2f', out_folder, rho);
mkdir(out_folder);
for m=1:2*nA,
    imwrite(outimg(:, :, :, m), sprintf('%s/%03d.jpg', out_folder, m),'jpg');
end