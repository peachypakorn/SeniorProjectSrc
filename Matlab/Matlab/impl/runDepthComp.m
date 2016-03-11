function runDepthComp(recompParam, fileID, startInd)
    global DATAPATH;
    global RESULTPATH;
    out_folder = sprintf('%s%s/depth', RESULTPATH, fileID.testName);
    mkdir(out_folder);

    maxDepth = fileID.maxDepth;

    if isstruct(recompParam),
        recompLowD = recompParam.recompLowD;
        recompMedi = recompParam.recompMedi;
        recompUpD = recompParam.recompUpD;
        recompImpD = recompParam.recompImpD;
        recompGtD = recompParam.recompGtD;
        recompHafD = recompParam.recompHafD;
        nocompHafD = recompParam.nocompHafD;
        compLowImpD = recompParam.compLowImpD;
    else,
        recompLowD = false;
        recompMedi = false;
        recompUpD = false;
        recompImpD = false;
        recompGtD = false;
        recompHafD = false;
        nocompHafD = false;
        compLowImpD = false;
    end

    tmpDepthP1 = 0;
    tmpDepthP2 = 0;
    tmpDepthPH1 = 0;
    tmpDepthPH2 = 0;
    tmpImgL1 = 0;
    tmpImgL2 = 0;
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
        nRn = floor(nR/2);
        nCn = floor(nC/2);
        nRnn = floor(nR/4);
        nCnn = floor(nC/4);
        imgL1 = imresize(img1, [nRn, nCn]);
        imgL2 = imresize(img2, [nRn, nCn]);
        imgLL1 = imresize(img1, [nRnn, nCnn]);
        imgLL2 = imresize(img2, [nRnn, nCnn]);

        lowD_file_L = sprintf('%s/lowD_%03d_L.png', out_folder, iImg);
        lowD_file_R = sprintf('%s/lowD_%03d_R.png', out_folder, iImg);
        if recompLowD || ~(exist(lowD_file_L, 'file') == 2 & exist(lowD_file_R, 'file') == 2),
            [depthLL1, depthLL2] = lagrDepthEstimate(imgLL1, imgLL2, ceil(maxDepth*0.25), 1e10);
            imwrite(0.5*(depthLL1/maxDepth + 1), lowD_file_L);
            imwrite(0.5*(depthLL2/maxDepth + 1), lowD_file_R);
        else,
            depthLL1 = round(double(imread(lowD_file_L)) * 2 * maxDepth / 255 - maxDepth);
            depthLL2 = round(double(imread(lowD_file_R)) * 2 * maxDepth / 255 - maxDepth);
        end

        medi_file_L = sprintf('%s/medi_%03d_L.png', out_folder, iImg);
        medi_file_R = sprintf('%s/medi_%03d_R.png', out_folder, iImg);
        if recompMedi || ~(exist(medi_file_L, 'file') == 2 & exist(medi_file_R, 'file') == 2),
            [depthLL1, depthLL2] = trilateralDepthFilter(depthLL1, depthLL2, imgLL1, imgLL2, 17, 25, 0.1, 2, 40);
            imwrite(0.5*(depthLL1/maxDepth + 1), medi_file_L);
            imwrite(0.5*(depthLL2/maxDepth + 1), medi_file_R);
        else,
            depthLL1 = round(double(imread(medi_file_L)) * 2 * maxDepth / 255 - maxDepth);
            depthLL2 = round(double(imread(medi_file_R)) * 2 * maxDepth / 255 - maxDepth);
        end

        upD_file_L = sprintf('%s/upD_%03d_L.png', out_folder, iImg);
        upD_file_R = sprintf('%s/upD_%03d_R.png', out_folder, iImg);
        if recompUpD || ~(exist(upD_file_L, 'file') == 2 & exist(upD_file_R, 'file') == 2),
            depthL1 = depthUpsample(depthLL1, 0, imgL1, imgLL1, 0, 5, 5, 0.1, 100);
            depthL2 = depthUpsample(depthLL2, 0, imgL2, imgLL2, 0, 5, 5, 0.1, 100);
            depthP1 = depthUpsample(depthL1, tmpDepthP1, img1, imgL1, tmpImgL1, 5, 5, 0.1, 100);
            depthP2 = depthUpsample(depthL2, tmpDepthP2, img2, imgL2, tmpImgL2, 5, 5, 0.1, 100);
            tmpDepthP1 = imresize(depthP1, [nRn, nCn]) * nCn / nC;
            tmpDepthP2 = imresize(depthP2, [nRn, nCn]) * nCn / nC;
            imwrite(0.5*(depthP1/maxDepth + 1), upD_file_L);
            imwrite(0.5*(depthP2/maxDepth + 1), upD_file_R);
        else,
            depthP1 = round(double(imread(upD_file_L)) * 2 * maxDepth / 255 - maxDepth);
            depthP2 = round(double(imread(upD_file_R)) * 2 * maxDepth / 255 - maxDepth);
        end

        impD_file_L = sprintf('%s/impD_%03d_L.png', out_folder, iImg);
        impD_file_R = sprintf('%s/impD_%03d_R.png', out_folder, iImg);
        if recompImpD || ~(exist(impD_file_L, 'file') == 2 & exist(impD_file_R, 'file') == 2),
            kw = ceil(sqrt(2*maxDepth)*6);
            depth1 = round(asymmetricSmooth(depthP1, 4, 4, maxDepth, maxDepth, kw, 6/maxDepth));
            depth2 = round(asymmetricSmooth(depthP2, 4, 4, maxDepth, maxDepth, kw, 6/maxDepth));
            imwrite(0.5*(depth1/maxDepth + 1), impD_file_L);
            imwrite(0.5*(depth2/maxDepth + 1), impD_file_R);
        else,
            depth1 = double(imread(impD_file_L)) * 2 * maxDepth / 255 - maxDepth;
            depth2 = double(imread(impD_file_R)) * 2 * maxDepth / 255 - maxDepth;
        end

        lowImpD_file_L = sprintf('%s/lowImpD_%03d_L.png', out_folder, iImg);
        lowImpD_file_R = sprintf('%s/lowImpD_%03d_R.png', out_folder, iImg);
        if iImg >= startInd && compLowImpD && ~(exist(lowImpD_file_L, 'file') == 2 & exist(lowImpD_file_R, 'file') == 2),
            kw = ceil(sqrt(2*maxDepth)*6);
            depthLI1 = round(asymmetricSmooth(imresize(depthLL1, [nR, nC]) * nC/nCnn, 4, 4, maxDepth, maxDepth, kw, 6/maxDepth));
            depthLI2 = round(asymmetricSmooth(imresize(depthLL2, [nR, nC]) * nC/nCnn, 4, 4, maxDepth, maxDepth, kw, 6/maxDepth));
            imwrite(0.5*(depthLI1/maxDepth + 1), lowImpD_file_L);
            imwrite(0.5*(depthLI2/maxDepth + 1), lowImpD_file_R);
        else,
        end

        gtOrig_file_L = sprintf('%s/gtOrig_%03d_L.png', out_folder, iImg);
        gtOrig_file_R = sprintf('%s/gtOrig_%03d_R.png', out_folder, iImg);
        gt_file_L = sprintf('%s/gtD_%03d_L.png', out_folder, iImg);
        gt_file_R = sprintf('%s/gtD_%03d_R.png', out_folder, iImg);
        gtImp_file_L = sprintf('%s/gtImpD_%03d_L.png', out_folder, iImg);
        gtImp_file_R = sprintf('%s/gtImpD_%03d_R.png', out_folder, iImg);
        if (exist(gtOrig_file_L, 'file') == 2) && (recompGtD || ...
                ~(exist(gt_file_L, 'file') == 2 && exist(gt_file_R, 'file') == 2)),
            depthGtOrig1 = imread(gtOrig_file_L);
            depthGt1 = round((double(depthGtOrig1) + fileID.gtAdd) * fileID.gtMul);
            depthGt1 = fillLRMax(depthGt1, fileID.gtAdd * fileID.gtMul);
            if (exist(gtOrig_file_R, 'file') == 2)
                depthGtOrig2 = imread(gtOrig_file_R);
                depthGt2 = round((double(depthGtOrig2) + fileID.gtAdd) * fileID.gtMul);
                depthGt2 = fillLRMax(depthGt2, fileID.gtAdd * fileID.gtMul);
            else,
                %depthGt2 = warpDisp(depthGt1, depthGtImp1, -1, 2, true);
            end

            %[depthF1, depthF2] = trilateralDepthFilter(depthGt1, depthGt2, img1, img2, 13, 9, 0.5, 1, 40);
            %depthGt1(depthGtOrig1 == 0) = depthF1(depthGtOrig1 == 0);
            depthGtImp1 = round(asymmetricSmooth(depthGt1, 0.25, 0.25, 10, 10, 11, 2));
            if (exist(gtOrig_file_R, 'file') == 2)
            %    depthGt2(depthGtOrig2 == 0) = depthF2(depthGtOrig2 == 0);
            else,
                depthGt2 = warpDisp(depthGt1, depthGtImp1, -1, 2, true);
            end
            depthGtImp2 = round(asymmetricSmooth(depthGt2, 0.25, 0.25, 10, 10, 11, 2));

            if cut < 0,
                depthGt1 = depthGt1(:, 1:nC+cut, :)-cut;
                depthGt1 = imresize(depthGt1, [nR, nC]) * (nC / (nC+cut));
                depthGtImp1 = depthGtImp1(:, 1:nC+cut, :)-cut;
                depthGtImp1 = imresize(depthGtImp1, [nR, nC]) * (nC / (nC+cut));
                depthGt2 = depthGt2(:, 1-cut:nC, :)-cut;
                depthGt2 = imresize(depthGt2, [nR, nC]) * (nC / (nC+cut));
                depthGtImp2 = depthGtImp2(:, 1-cut:nC, :)-cut;
                depthGtImp2 = imresize(depthGtImp2, [nR, nC]) * (nC / (nC+cut));
            elseif cut > 0,
                depthGt1 = depthGt1(:, 1+cut:nC, :)-cut;
                depthGt1 = imresize(depthGt1, [nR, nC]) * (nC / (nC-cut));
                depthGtImp1 = depthGtImp1(:, 1+cut:nC, :)-cut;
                depthGtImp1 = imresize(depthGtImp1, [nR, nC]) * (nC / (nC-cut));
                depthGt2 = depthGt2(:, 1:nC-cut, :)-cut;
                depthGt2 = imresize(depthGt2, [nR, nC]) * (nC / (nC-cut));
                depthGtImp2 = depthGtImp2(:, 1:nC-cut, :)-cut;
                depthGtImp2 = imresize(depthGtImp2, [nR, nC]) * (nC / (nC-cut));
            end

            imwrite(0.5*(depthGt1/maxDepth + 1), gt_file_L);
            imwrite(0.5*(depthGt2/maxDepth + 1), gt_file_R);
            imwrite(0.5*(depthGtImp1/maxDepth + 1), gtImp_file_L);
            imwrite(0.5*(depthGtImp2/maxDepth + 1), gtImp_file_R);
        end

        if ~nocompHafD,
            hafD_file_L = sprintf('%s/hafD_%03d_L.png', out_folder, iImg);
            hafD_file_R = sprintf('%s/hafD_%03d_R.png', out_folder, iImg);
            if recompHafD || ~(exist(hafD_file_L, 'file') == 2 & exist(hafD_file_R, 'file') == 2),
                [depthL1, depthL2] = lagrDepthEstimate(imgL1, imgL2, ceil(maxDepth*0.5), 1e10);
                imwrite(0.5*(depthL1/maxDepth + 1), hafD_file_L);
                imwrite(0.5*(depthL2/maxDepth + 1), hafD_file_R);
            else,
                depthL1 = round(double(imread(hafD_file_L)) * 2 * maxDepth / 255 - maxDepth);
                depthL2 = round(double(imread(hafD_file_R)) * 2 * maxDepth / 255 - maxDepth);
            end

            medi_file_L = sprintf('%s/mediH_%03d_L.png', out_folder, iImg);
            medi_file_R = sprintf('%s/mediH_%03d_R.png', out_folder, iImg);
            if recompMedi || recompHafD || ~(exist(medi_file_L, 'file') == 2 & exist(medi_file_R, 'file') == 2),
                [depthL1, depthL2] = trilateralDepthFilter(depthL1, depthL2, imgL1, imgL2, 17, 25, 0.1, 2, 20);
                imwrite(0.5*(depthL1/maxDepth + 1), medi_file_L);
                imwrite(0.5*(depthL2/maxDepth + 1), medi_file_R);
            else,
                depthL1 = round(double(imread(medi_file_L)) * 2 * maxDepth / 255 - maxDepth);
                depthL2 = round(double(imread(medi_file_R)) * 2 * maxDepth / 255 - maxDepth);
            end

            upD_file_L = sprintf('%s/upDH_%03d_L.png', out_folder, iImg);
            upD_file_R = sprintf('%s/upDH_%03d_R.png', out_folder, iImg);
            if recompUpD || recompHafD || ~(exist(upD_file_L, 'file') == 2 & exist(upD_file_R, 'file') == 2),
                depthP1 = depthUpsample(depthL1, tmpDepthPH1, img1, imgL1, tmpImgL1, 5, 5, 0.1, 100);
                depthP2 = depthUpsample(depthL2, tmpDepthPH2, img2, imgL2, tmpImgL2, 5, 5, 0.1, 100);
                tmpDepthPH1 = imresize(depthP1, [nRn, nCn]) * nCn / nC;
                tmpDepthPH2 = imresize(depthP2, [nRn, nCn]) * nCn / nC;
                imwrite(0.5*(depthP1/maxDepth + 1), upD_file_L);
                imwrite(0.5*(depthP2/maxDepth + 1), upD_file_R);
            else,
                depthP1 = round(double(imread(upD_file_L)) * 2 * maxDepth / 255 - maxDepth);
                depthP2 = round(double(imread(upD_file_R)) * 2 * maxDepth / 255 - maxDepth);
            end

            impD_file_L = sprintf('%s/impDH_%03d_L.png', out_folder, iImg);
            impD_file_R = sprintf('%s/impDH_%03d_R.png', out_folder, iImg);
            if recompImpD || recompHafD || ~(exist(impD_file_L, 'file') == 2 & exist(impD_file_R, 'file') == 2),
                kw = ceil(sqrt(2*maxDepth)*6);
                depth1 = round(asymmetricSmooth(depthP1, 4, 4, maxDepth, maxDepth, kw, 6/maxDepth));
                depth2 = round(asymmetricSmooth(depthP2, 4, 4, maxDepth, maxDepth, kw, 6/maxDepth));
                imwrite(0.5*(depth1/maxDepth + 1), impD_file_L);
                imwrite(0.5*(depth2/maxDepth + 1), impD_file_R);
            else,
                depth1 = double(imread(impD_file_L)) * 2 * maxDepth / 255 - maxDepth;
                depth2 = double(imread(impD_file_R)) * 2 * maxDepth / 255 - maxDepth;
            end
        end
        tmpImgL1 = imgL1;
        tmpImgL2 = imgL2;
    end
end