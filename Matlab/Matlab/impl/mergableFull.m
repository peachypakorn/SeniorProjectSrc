
%source_image_1 = imresize(double(imread('data/bbb_106_left.jpg')) / 255, 0.5);
%source_image_2 = imresize(double(imread('data/bbb_106_right.jpg')) / 255, 0.5);

%source_image_1 = imresize(double(imread('data/bbb_252p7_left.jpg')) / 255, 0.5);
%source_image_2 = imresize(double(imread('data/bbb_252p7_right.jpg')) / 255, 0.5);

%source_image_1 = repmat(double(imread('data/check_1.png')) / 255, [1, 1, 3]);
%source_image_2 = repmat(double(imread('data/check_2.png')) / 255, [1, 1, 3]);

%source_image_1 = double(imread('data/check7_1.png')) / 255;
%source_image_2 = double(imread('data/check7_2.png')) / 255;

%source_image_1 = double(imread('data/window1_1.jpg')) / 255;
%source_image_2 = double(imread('data/window1_2.jpg')) / 255;

%source_image_1 = double(imread('data/gift_1.jpg')) / 255;
%source_image_2 = double(imread('data/gift_2.jpg')) / 255;

%source = imresize(double(imread('data/bicycle.jpg')) / 255, 1/3);
%source_image_1 = source(:, 1:size(source, 2)/2, :);
%source_image_2 = source(:, size(source, 2)/2+1:size(source, 2), :);

source_image_1 = double(imread('data/frame1.jpg')) / 255;
source_image_2 = double(imread('data/frame2.jpg')) / 255;
%source_image_1 = ones([32, 128, 3]);
%source_image_2 = ones([32, 128, 3]);

%source_image_1 = imresize(double(imread('data/art_1.png')) / 255, 0.5) ;
%source_image_1 = source_image_1(:, 1:1024, :);
%source_image_2 = imresize(double(imread('data/art_2.png')) / 255, 0.5) ;
%source_image_2 = source_image_2(:, 1:1024, :);
%depth_image = (double(imread('data/art_d.png'))) / 4;
%depth_image = depth_image(:, 1:1024, 1);

%source_image_1 = (double(imread('data/doll_1.png')) / 255) ;
%source_image_1 = source_image_1(:, 1:512, :);
%source_image_2 = (double(imread('data/doll_2.png')) / 255) ;
%source_image_2 = source_image_2(:, 1:512, :);

%source_image_1 = double(imread('data/transparency_left.jpg')) / 255;
%source_image_2 = double(imread('data/transparency_right.jpg')) / 255;

%source_image_1 = imresize(double(imread('data/camera_left.jpg')) / 255, 0.85);
%source_image_2 = imresize(double(imread('data/camera_right.jpg')) / 255, 0.85);
%source_image_1 = source_image_1(1:128, :, :);
%source_image_2 = source_image_2(1:128, :, :);

nC = size(source_image_1, 2);
po2 = 2^ceil(log2(nC));
nCL = floor((po2 - nC)/2);
nCR = po2 - nC - nCL;

source_image_1 = [source_image_1(:, nCL:-1:1, :), source_image_1(:, :, :), source_image_1(:, nC:-1:nC-nCR+1, :)];
source_image_2 = [source_image_2(:, nCL:-1:1, :), source_image_2(:, :, :), source_image_2(:, nC:-1:nC-nCR+1, :)];

amp = [0 1 2 3];
aamp = [-4 -3 -2 -1 0 1 2 3];
%amp = ((1:2:23) - 8) / 8;
%aamp = ((-23:2:23) - 8) / 8;
%amp = ((1:2:23)) / 6 - 1;
%aamp = ((-23:2:23))  / 6 - 1;
nA = length(amp);
nAA = length(aamp);
trim = 8;
mmtrim = 0;
rho = 0.00;
rhor = 0.0;
outimg = zeros([ size(source_image_1) 2*nA ]);
%debug only
debugStrip = 128;
depthFactor = 32;
depthOffset = 0.5;

[outimg(:, :, :, nAA:-1:1), ldepth, lmovPyr] = mergableMultiImageSynthDiamond(...
        source_image_1(:, :, :), source_image_2(:, :, :), aamp, debugStrip, trim, depthFactor, depthOffset, rho, rhor, mmtrim, 0, 0);
[tmpimg, rdepth, rmovPyr] = mergableMultiImageSynthDiamond(...
        source_image_2(:, :, :), source_image_1(:, :, :), aamp, debugStrip, trim, depthFactor, 1-depthOffset, rho, rhor, mmtrim, 0, 0);

outimg(:, :, :, nAA+1:2*nA) = tmpimg(:, :, :, 2*nAA-2*nA+1:nAA);
for aid = 1:nAA-nA,
    v1 = abs(-aamp(aid));
    v2 = abs(1+aamp(aid));
    outimg(:, :, :, nA+aid) = outimg(:, :, :, nA+aid) * (v2 / (v1+v2)) + tmpimg(:, :, :, nAA-nA+aid) * (v1 / (v1+v2));
    outimg(:, :, :, nA+1-aid) = outimg(:, :, :, nA+1-aid) * (v1 / (v1+v2)) + tmpimg(:, :, :, nAA-nA+1-aid) * (v2 / (v1+v2));
end

out_folder = sprintf('result/syn_%.2f', rho);
if rhor > 0,
    out_folder = sprintf('%s_%.2f', out_folder, rhor);
end
    
mkdir(out_folder);
for m=1:2*nA,
    imwrite(outimg(:, nCL+1:nCL+nC, :, m), sprintf('%s/%03d.png', out_folder, m),'png');
end
imwrite(ldepth / depthFactor + depthOffset, sprintf('%s/depth_l.png', out_folder),'png');
imwrite(-rdepth / depthFactor + 1-depthOffset, sprintf('%s/depth_r.png', out_folder),'png');
