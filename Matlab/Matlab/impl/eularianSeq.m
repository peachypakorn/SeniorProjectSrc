nIvl = 2;
%in_seq = 'data/bbb_seq_305/bbb_305';
in_seq = 'data/dracular-seq-215/ud-';
%in_seq = 'data/airport-seq-243/airport-243-ud-pa';
nSeq = 24;
nStart = 215;
mult = 2/3;
maxDLagr = 12;
dsrLagr = 0.25;

rho = 0.0;
ht = 6;
%out_folder = sprintf('result/eularian-1D/bbb-seq-305/forcematch_eulaOlagr_v2_%.2f', rho);
out_folder = sprintf('result/eularian-1D/dracular-seq-215/eulaOlagr_v2_%.2f', rho);
%out_folder = sprintf('result/lagrangian/dracular-seq-215/directwarp_%.2f', rho);
mkdir(out_folder);

for seqI = (nStart-1) + (1:nIvl:nSeq),
    seqI
    t0 = cputime();
    %source_image_1 = imresize(double(imread(sprintf('%s_right_%03d.jpg', in_seq, seqI))) / 255, mult);
    %source_image_2 = imresize(double(imread(sprintf('%s_left_%03d.jpg', in_seq, seqI))) / 255, mult);
    source_image_1 = imresize(double(imread(sprintf('%s%04d-0.jpg', in_seq, seqI))) / 255, mult);
    source_image_2 = imresize(double(imread(sprintf('%s%04d-1.jpg', in_seq, seqI))) / 255, mult);
    %source_image = imresize(double(imread(sprintf('%s_%03d.jpg', in_seq, seqI))) / 255, mult);
    %nC = size(source_image, 2);
    %source_image_1 = source_image(:, 1:nC/2, :);
    %source_image_2 = source_image(:, nC/2+1:nC, :);

    %outimg = eularian1D(source_image_1, source_image_2, rho, ht, 4, 1, 3, 0);
    outimg = eulaOnLagr(source_image_1, source_image_2, rho, ht, 4, 1, dsrLagr, maxDLagr*dsrLagr, 20, 0, 0, 0);
    %[outimg, depthL, depthR] = lagrangian(source_image_1, source_image_2, 4, 1, 1, maxDLagr, 80, 0, 0);

    for m=1:8
        imwrite(outimg(:, :, :, m), sprintf('%s/%03d_%03d.jpg', out_folder, seqI, m),'jpg');
    end
    %imwrite(0.5*depthL/maxDLagr + 0.5, sprintf('%s/%03d_depthL.jpg', out_folder, seqI),'jpg');
    %imwrite(0.5*depthR/maxDLagr + 0.5, sprintf('%s/%03d_depthR.jpg', out_folder, seqI),'jpg');
    cputime() - t0
end