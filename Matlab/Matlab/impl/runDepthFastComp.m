function runDepthFastComp(recompParam, fileID, startInd)
    global DATAPATH;
    global RESULTPATH;
    out_folder = sprintf('%s%s/depth', RESULTPATH, fileID.testName);
    mkdir(out_folder);

    if isstruct(recompParam),
        recompLowFastD = recompParam.recompLowFastD;
    else,
        recompLowFastD = false;
    end

    maxDepth = fileID.maxDepth;
    prevVals = struct();

    for iImg = 1:fileID.nImg,
        disp(sprintf('%s-depth-%d', fileID.testName, iImg));
        t0 = cputime();
        if strcmp(fileID.inpType, 'separate'),
            img1 = double(imread(sprintf('%s%s%s', DATAPATH, fileID.inpPrefix, fileID.inpImg1{iImg}))) / 255;
            img2 = double(imread(sprintf('%s%s%s', DATAPATH, fileID.inpPrefix, fileID.inpImg2{iImg}))) / 255;
        elseif strcmp(fileID.inpType, 'combined'),
            img = double(imread(sprintf('%s%s%s', DATAPATH, fileID.inpPrefix, fileID.inpImg{iImg}))) / 255;
            nC = size(img, 2);
            img1 = img(:, 1:nC/2, :);
            img2 = img(:, nC/2+1:nC, :);
        elseif strcmp(fileID.inpType, 'combinedInv'),
            img = double(imread(sprintf('%s%s%s', DATAPATH, fileID.inpPrefix, fileID.inpImg{iImg}))) / 255;
            nC = size(img, 2);
            img2 = img(:, 1:nC/2, :);
            img1 = img(:, nC/2+1:nC, :);
        elseif strcmp(fileID.inpType, 'topDown'),
            img = double(imread(sprintf('%s%s%s', DATAPATH, fileID.inpPrefix, fileID.inpImg{iImg}))) / 255;
            nR = size(img, 1);
            img1 = img(1:nR/2, :, :);
            img2 = img(nR/2+1:nR, :, :);
        else,
            error('unknown inpType');
        end

        if isfield(fileID, 'resize'),
            img1 = imresize(img1, fileID.resize);
            img2 = imresize(img2, fileID.resize);
        end
        if isfield(fileID, 'cut'),
            cut = fileID.cut;
            nC = size(img1, 2);
            if cut < 0,
                img1 = img1(:, 1:nC+cut, :);
                img2 = img2(:, 1-cut:nC, :);
            elseif cut > 0,
                img1 = img1(:, 1+cut:nC, :);
                img2 = img2(:, 1:nC-cut, :);
            end
            img1 = imresize(img1, [size(img1, 1), nC]);
            img2 = imresize(img2, [size(img1, 1), nC]);
        else,
            cut = 0;
        end

        nR = size(img1, 1);
        nC = size(img1, 2);

        lowFastD_file_L = sprintf('%s/lowFastD_%03d_L.png', out_folder, iImg);
        lowFastD_file_R = sprintf('%s/lowFastD_%03d_R.png', out_folder, iImg);
        if recompLowFastD || ~(exist(lowFastD_file_L, 'file') == 2 & exist(lowFastD_file_R, 'file') == 2),
            kR = ceil(nR / 350) * 2 + 1;
            kC = ceil(nC / 300) * 2 + 1;
            [depthLL1, depthLL2, tmp1, tmp2, prevVals] = lagrDepthFast(img1, img2, -maxDepth, 2, maxDepth+1, 4, 4, kR, kC, kR, kC, 1, prevVals, 0.5);
            depthLL1 = maxDepth - 2*depthLL1;
            depthLL2 = maxDepth - 2*depthLL2;

            imwrite(0.5*(1+depthLL1/maxDepth), lowFastD_file_L);
            imwrite(0.5*(1+depthLL2/maxDepth), lowFastD_file_R);
        else,
            depthLL1 = double(imread(lowFastD_file_L)) * 2 * maxDepth / 255 - maxDepth;
            depthLL2 = double(imread(lowFastD_file_R)) * 2 * maxDepth / 255 - maxDepth;
        end
    end
end