nIvl = 1;
nSeq = 23;
nStart = 25;
in_seq = 'data/sintel-seq-290/sintel-290';
out_seq = 'data/sintel-seq-290/sintel-290-ud';

htor = 20;
vtor = 5;
crop = 2;
mult = 1;

for seqI = (nStart-1) + (1:nIvl:nSeq),
    seqI
    source_image = imresize(double(imread(sprintf('%s_%03d.jpg', in_seq, seqI))) / 255, mult);
    source_image_p = imresize(double(imread(sprintf('%s_%03d.jpg', in_seq, seqI+1))) / 255, mult);
    nC = size(source_image, 2);
    source_image_2 = source_image_p(:, 1:nC/2, :);
    source_image_1 = source_image(:, nC/2+1:nC, :);

    nR = size(source_image_1, 1);
    nC = size(source_image_1, 2);

    [imgC1, imgC2] = verticalDistortion(source_image_1, source_image_2, htor, vtor);
    imgC1 = imresize(imgC1(crop+1:nR-crop, :, :), [nR, nC]);
    imgC2 = imresize(imgC2(crop+1:nR-crop, :, :), [nR, nC]);

    imgC = zeros(nR, nC*2, size(source_image_1, 3));
    imgC(:, 1:nC, :) = imgC1;
    imgC(:, nC+1:2*nC, :) = imgC2;
    imwrite(imgC, sprintf('%s_%03d.jpg', out_seq, seqI), 'jpg');
end
