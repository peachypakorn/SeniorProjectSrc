
%source_image_1 = double(imread('../bbb_106_left.jpg')) / 255;
%source_image_2 = double(imread('../bbb_106_right.jpg')) / 255;
%source_image_1 = source_image_1(:, 801:1824, :);
%source_image_2 = source_image_2(:, 801:1824, :);


%source_image_1 = double(imread('../bbb_252p7_left.jpg')) / 255;
%source_image_2 = double(imread('../bbb_252p7_right.jpg')) / 255;
%source_image_1 = source_image_1(:, 301:1324, :);
%source_image_2 = source_image_2(:, 301:1324, :);

source_image_1 = double(imread('../frame1.jpg')) / 255;
source_image_2 = double(imread('../frame2.jpg')) / 255;
nC = size(source_image_1, 2);
source_image_1 = [source_image_1(:, :, :), zeros([size(source_image_1, 1), 1024 - nC, 3])];
source_image_2 = [source_image_2(:, :, :), zeros([size(source_image_1, 1), 1024 - nC, 3])];

%source_image_1 = (double(imread('../art_1.png')) / 255) ;
%source_image_1 = source_image_1(:, 1:1024, :);
%source_image_2 = (double(imread('../art_2.png')) / 255) ;
%source_image_2 = source_image_2(:, 1:1024, :);
%depth_image = (double(imread('../art_d.png'))) / 4;
%depth_image = depth_image(:, 1:1024, 1);

%source_image_1 = (double(imread('../doll_1.png')) / 255) ;
%source_image_1 = source_image_1(:, 1:512, :);
%source_image_2 = (double(imread('../doll_2.png')) / 255) ;
%source_image_2 = source_image_2(:, 1:512, :);

amp = [-1 -2 -3];
oc = 16;
outimg = zeros([ size(source_image_1) length(amp) ]);

%out_folder = sprintf('tmp_debug');
%mkdir(out_folder);
%imwrite(repmat(source_image_2, [1, 1, 1]), sprintf('%s/001.png', out_folder),'png');
%imwrite(repmat(source_image_1, [1, 1, 1]), sprintf('%s/002.png', out_folder),'png');

%outimg = mergableMultiImageSynth(source_image_1(:, :, :), source_image_2(:, :, :), amp, oc, 128, 32, 32, 64, 0.0);
outimg = mergableMultiImageSynthDiamond(source_image_1(:, :, :), source_image_2(:, :, :), amp, 32, 8, 32, 0.5, 0.0);
%outimg = mergableMultiImageSynthFromDepth(source_image_1(:, :, :), depth_image, amp, oc, 32, 8);

out_folder = sprintf('syn_mergable_bilateral');
mkdir(out_folder);
imwrite(repmat(source_image_2(:, 1:nC, :), [1, 1, 1]), sprintf('%s/001.png', out_folder),'png');
imwrite(repmat(source_image_1(:, 1:nC, :), [1, 1, 1]), sprintf('%s/002.png', out_folder),'png');
for m=1:length(amp),
    %imwrite(outimg(:, :, :, m), sprintf('%s/%03d.png', out_folder, m),'png');
    imwrite(repmat(outimg(:, 1:nC, :, m), [1, 1, 1]), sprintf('%s/%03d.png', out_folder, m+2),'png');
end
