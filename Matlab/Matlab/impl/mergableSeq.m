nIvl = 2;

amp = [-1 1 3 5] / 4.0;
aamp = [-9 -7 -5 -3 -1 1 3 5] / 4.0;
nA = length(amp);
nAA = length(aamp);
mmtrim = 0;
trim = 16;
rho = 0.5;
rhor = 0.2;
prvAlpha = 0.5;

in_seq = 'data/bbb_seq_69/bbb_69';
out_folder = sprintf('result/bbb_seq_69_tsmooth2_reduce_%.2f', rho);
nSeq = 143;
mult = 2.0/3;

%in_seq = 'data/bbb_seq_155/bbb_155';
%out_folder = sprintf('result/bbb_seq_155_tsmooth2_reduce_%.2f', rho);
%nSeq = 87;
%mult = 1.0/3;

if rhor > 0,
    out_folder = sprintf('%s_%.2f', out_folder, rhor);
end

mkdir(out_folder);

%debug only
debugStrip = 0;
depthFactor = 64;
depthOffset = 0.5;

lmovPyr = 0;
rmovPyr = 0;

for seqI = 111:nIvl:nSeq,
    seqI
    source_image_1 = imresize(double(imread(sprintf('%s_left_%03d.jpg', in_seq, seqI))) / 255, mult);
    source_image_2 = imresize(double(imread(sprintf('%s_right_%03d.jpg', in_seq, seqI))) / 255, mult);


    nC = size(source_image_1, 2);
    po2 = 2^ceil(log2(nC));
    nCL = floor((po2 - nC)/2);
    nCR = po2 - nC - nCL;

    source_image_1 = [source_image_1(:, nCL:-1:1, :), source_image_1(:, :, :), source_image_1(:, nC:-1:nC-nCR+1, :)];
    source_image_2 = [source_image_2(:, nCL:-1:1, :), source_image_2(:, :, :), source_image_2(:, nC:-1:nC-nCR+1, :)];

    outimg = zeros([ size(source_image_1) 2*nA ]);

    [outimg(:, :, :, nAA:-1:1), ldepth, lmovPyr] = mergableMultiImageSynthDiamond(...
            source_image_1(:, :, :), source_image_2(:, :, :), aamp, debugStrip, trim, depthFactor, depthOffset, rho, rhor, mmtrim, lmovPyr, prvAlpha);
    for m=1:nA+nA-nAA,
        imwrite(outimg(:, nCL+1:nCL+nC, :, m), sprintf('%s/%03d_%03d.jpg', out_folder, seqI, m),'jpg');
    end
    imwrite(ldepth / depthFactor + depthOffset, sprintf('%s/%03d_dl.jpg', out_folder, seqI),'jpg');

    [tmpImgs, rdepth, rmovPyr] = mergableMultiImageSynthDiamond(...
            source_image_2(:, :, :), source_image_1(:, :, :), aamp, debugStrip, trim, depthFactor, 1-depthOffset, rho, rhor, mmtrim, rmovPyr, prvAlpha);
    outimg(:, :, :, nAA+1:2*nA) = tmpImgs(:, :, :, 2*nAA-2*nA+1:nAA);
    for aid = 1:nAA-nA,
        v1 = abs(-aamp(aid));
        v2 = abs(1+aamp(aid));
        outimg(:, :, :, nA+aid) = outimg(:, :, :, nA+aid) * (v2 / (v1+v2)) + tmpImgs(:, :, :, nAA-nA+aid) * (v1 / (v1+v2));
        outimg(:, :, :, nA+1-aid) = outimg(:, :, :, nA+1-aid) * (v1 / (v1+v2)) + tmpImgs(:, :, :, nAA-nA+1-aid) * (v2 / (v1+v2));
    end
    for m=nA+nA-nAA+1:2*nA,
        imwrite(outimg(:, nCL+1:nCL+nC, :, m), sprintf('%s/%03d_%03d.jpg', out_folder, seqI, m),'jpg');
    end
    imwrite(-rdepth / depthFactor + 1-depthOffset, sprintf('%s/%03d_dr.jpg', out_folder, seqI),'jpg');
end
