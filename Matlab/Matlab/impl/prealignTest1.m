nIvl = 1;
nSeq = 24;
nStart = 721;
%in_seq = 'data/heidelberg-seq-338/heidelberg-338-ud';
%out_seq = 'data/heidelberg-seq-338/heidelberg-338-ud-pa';
in_seq = 'data/demo-seq-721/';
out_seq = 'data/demo-seq-721/pa-';

maxD = 20;
preD = -maxD-1;
alpha = 0.8;

for seqI = (nStart-1) + (1:nIvl:nSeq),
    seqI
    %source_image = double(imread(sprintf('%s_%03d.jpg', in_seq, seqI))) / 255;
    %nC = size(source_image, 2);
    %source_image_1 = source_image(:, 1:nC/2, :);
    %source_image_2 = source_image(:, nC/2+1:nC, :);
    source_image_1 = double(imread(sprintf('%s%04d-0.jpg', in_seq, seqI))) / 255;
    source_image_2 = double(imread(sprintf('%s%04d-1.jpg', in_seq, seqI))) / 255;

    nR = size(source_image_1, 1);
    nC = size(source_image_1, 2);

    [imgC1, imgC2, preD] = prealign(source_image_1, source_image_2, maxD, preD, alpha);

    nCC = nC - maxD;
    imgC = zeros(nR, nCC*2, size(source_image_1, 3));
    imgC(:, 1:nCC, :) = imgC1;
    imgC(:, nCC+1:2*nCC, :) = imgC2;
    %imwrite(imgC, sprintf('%s_%03d.jpg', out_seq, seqI), 'jpg');
    imwrite(imgC1, sprintf('%s%04d-0.jpg', out_seq, seqI), 'jpg');
    imwrite(imgC2, sprintf('%s%04d-1.jpg', out_seq, seqI), 'jpg');
end
