
%source_image_1 = double(imread('data/frame1.jpg')) / 255;
%source_image_2 = double(imread('data/frame2.jpg')) / 255;
%source_image_1 = double(imread('data/bbb_69_left_062.jpg')) / 255;
%source_image_2 = double(imread('data/bbb_69_right_062.jpg')) / 255;
%source_image_1 = double(imread('data/bbb_252p7_right.jpg')) / 255;
%source_image_2 = double(imread('data/bbb_252p7_left.jpg')) / 255;
%source_image_1 = double(imread('data/bbb_305_right_228.jpg')) / 255;
%source_image_2 = double(imread('data/bbb_305_left_228.jpg')) / 255;
%source_image = double(imread('data/video_2_aligned_1.jpg')) / 255;
%source_image_1 = source_image(:, size(source_image, 2)/2+1:size(source_image, 2), :);
%source_image_2 = source_image(:, 1:size(source_image, 2)/2, :);

source_image_1 = double(imread('data/transparency_right.jpg')) / 255;
source_image_2 = double(imread('data/transparency_left.jpg')) / 255;
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

htor = 30;
vtor = 10;
crop = 10;

nR = size(source_image_1, 1);
nC = size(source_image_1, 2);

[imgC1, imgC2] = verticalDistortion(source_image_1, source_image_2, htor, vtor);
imgC1 = imresize(imgC1(crop+1:nR-crop, :, :), [nR, nC]);
imgC2 = imresize(imgC2(crop+1:nR-crop, :, :), [nR, nC]);
imgC = zeros(nR, nC*2, size(source_image_1, 3));
imgC(:, 1:nC, :) = imgC1;
imgC(:, nC+1:2*nC, :) = imgC2;
imwrite(imgC, 'data/transparency_undist.jpg', 'jpg');
