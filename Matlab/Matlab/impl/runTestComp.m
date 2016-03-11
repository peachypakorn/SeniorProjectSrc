function runTestComp(fileID, testParam, runParam, methodID, startInd, dsId, depthRecompParam, noRecomp)
    global RESULTPATH;
    global DATAPATH;
    out_folder = sprintf('%s%s/%s/%s-%s', RESULTPATH, fileID.testName, ...
            runParam.runName, getMethodNameComp(methodID), getDepthSelect(dsId));
    mkdir(out_folder);

    phsCmpO = 0;
    for iImg = startInd:fileID.nImg,
        disp(sprintf('%s-%s-%s-%d-%s', fileID.testName, runParam.runName, getMethodNameComp(methodID), iImg, getDepthSelect(dsId)));
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

        close all;

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
        end
        depthImg1 = 0;
        depthImg2 = 0;
        depth_folder = sprintf('%s%s/depth', RESULTPATH, fileID.testName);
        if dsId ~= 0,
            depth_file1 = sprintf('%s/%s_%03d_L.png', depth_folder, getDepthSelect(dsId), iImg);
            depth_file2 = sprintf('%s/%s_%03d_R.png', depth_folder, getDepthSelect(dsId), iImg);
            if dsId < 9,
                if ~(exist(depth_file1, 'file') == 2 & exist(depth_file2, 'file') == 2),
                    disp('no depth image found');
                    runDepthComp(depthRecompParam, fileID, iImg);
                elseif isstruct(depthRecompParam) & iImg == startInd,
                    disp('recomputing depth');
                    runDepthComp(depthRecompParam, fileID, iImg);
                end
            else,
                if ~(exist(depth_file1, 'file') == 2 & exist(depth_file2, 'file') == 2),
                    disp('no depth image found');
                    runDepthFastComp(depthRecompParam, fileID, iImg);
                elseif isstruct(depthRecompParam) & iImg == startInd,
                    disp('recomputing depth');
                    runDepthFastComp(depthRecompParam, fileID, iImg);
                end
            end
            if exist(depth_file1, 'file') == 2 & exist(depth_file2, 'file') == 2,
                depthImg1 = double(imread(depth_file1)) * 2 * fileID.maxDepth / 255 - fileID.maxDepth;
                depthImg2 = double(imread(depth_file2)) * 2 * fileID.maxDepth / 255 - fileID.maxDepth;
            else
                error('depth map error');
            end
        end

        if iImg == startInd && ((isfield(runParam, 'lincomp') && runParam.lincomp) || ...
            (isfield(runParam, 'remap') && runParam.remap)),
            depth_file1 = sprintf('%s/lowD_%03d_L.png', depth_folder, 1);
            depthImg = double(imread(depth_file1)) * 2 * fileID.maxDepth / 255 - fileID.maxDepth;
            minD = prctile(depthImg(:), 1);
            maxD = prctile(depthImg(:), 99);
            fact = size(img1, 2) / size(depthImg, 2);
            if (isfield(runParam, 'lincomp') && runParam.lincomp),
                runParam.dA = (40 / (fact * (maxD-minD))) / (2*runParam.nA-1);
                cut = round(round((minD + maxD)/2) * fact);
            else,
                cut = 0;
            end
            if (isfield(runParam, 'remap') && runParam.remap),
                runParam.origMinD = min(0, minD * fact - cut);
                runParam.origMaxD = max(0, maxD * fact - cut);
            end
        else,
            cut = 0;
        end

        depthImgNonFilt1 = 0;
        depthImgNonFilt2 = 0;
        if strcmp(getDepthSelect(dsId), 'impD'),
            depth_file1 = sprintf('%s/%s_%03d_L.png', depth_folder, 'upD', iImg);
            depth_file2 = sprintf('%s/%s_%03d_R.png', depth_folder, 'upD', iImg);
            if exist(depth_file1, 'file') == 2 & exist(depth_file2, 'file') == 2,
                depthImgNonFilt1 = double(imread(depth_file1)) * 2 * fileID.maxDepth / 255 - fileID.maxDepth;
                depthImgNonFilt2 = double(imread(depth_file2)) * 2 * fileID.maxDepth / 255 - fileID.maxDepth;
            else
                error('non-filtered depth map error');
            end
        elseif strcmp(getDepthSelect(dsId), 'gtImpD'),
            depth_file1 = sprintf('%s/%s_%03d_L.png', depth_folder, 'gtD', iImg);
            depth_file2 = sprintf('%s/%s_%03d_R.png', depth_folder, 'gtD', iImg);
            if exist(depth_file1, 'file') == 2 & exist(depth_file2, 'file') == 2,
                depthImgNonFilt1 = double(imread(depth_file1)) * 2 * fileID.maxDepth / 255 - fileID.maxDepth;
                depthImgNonFilt2 = double(imread(depth_file2)) * 2 * fileID.maxDepth / 255 - fileID.maxDepth;
            else
                error('non-filtered depth map error');
            end
        else,
            depthImgNonFilt1 = depthImg1;
            depthImgNonFilt2 = depthImg2;
        end

        nCOrig = size(img1, 2);
        if cut < 0,
            nC = size(img1, 2);
            nCd = size(depthImg1, 2);
            cutD = round(cut * nCd / nC);
            img1 = img1(:, 1:nC+cut, :);
            depthImg1 = depthImg1(:, 1:nCd+cutD, :)-cutD;
            depthImgNonFilt1 = depthImgNonFilt1(:, 1:nCd+cutD, :)-cutD;
            img2 = img2(:, 1-cut:nC, :);
            depthImg2 = depthImg2(:, 1-cutD:nCd, :)-cutD;
            depthImgNonFilt2 = depthImgNonFilt2(:, 1-cutD:nCd, :)-cutD;
        elseif cut > 0,
            nC = size(img1, 2);
            nCd = size(depthImg1, 2);
            cutD = round(cut * nCd / nC);
            img1 = img1(:, 1+cut:nC, :);
            depthImg1 = depthImg1(:, 1+cutD:nCd, :)-cutD;
            depthImgNonFilt1 = depthImgNonFilt1(:, 1+cutD:nCd, :)-cutD;
            img2 = img2(:, 1:nC-cut, :);
            depthImg2 = depthImg2(:, 1:nCd-cutD, :)-cutD;
            depthImgNonFilt2 = depthImgNonFilt2(:, 1:nCd-cutD, :)-cutD;
        end


        if methodID == 1 || strcmp(methodID, 'eol'),
            mem = 0;
            nA = runParam.nA;
            nR = size(img1, 1);
            nC = size(img1, 2);

            for vId = 1:2*nA,
                disp(sprintf('View %d: ', vId));
                out_file = sprintf('%s/%03d_%03d.jpg', out_folder, iImg, vId);
                if noRecomp && exist(out_file, 'file') == 2,
                    disp('not recomputing');
                    continue;
                end
                [synth, mem, phsCmpO] = eulaOnLagrComp(reflectEx(img1), reflectEx(img2), reflectEx(depthImg1), reflectEx(depthImg2), testParam, runParam, mem, vId, false, phsCmpO);
                nCRef = size(synth, 2);
                nEx = (nCRef - nC)/2;
                outimg = synth(:, nEx+(1:nC), :);
                imwrite(imresize(outimg, [nR, nCOrig]), out_file,'jpg');
            end
        elseif methodID == 2 || strcmp(methodID, 'lagr') || methodID == 7 || strcmp(methodID, 'lagrSmooth'),
            nA = runParam.nA;
            dA = runParam.dA;

            dspLagr = size(depthImg1, 2) / size(img1, 2);
            nR = size(img1, 1);
            nC = size(img1, 2);
            nCol = size(img1, 3);
            rho = runParam.rho;
            depthImg1 = imresize(depthImg1, [nR, nC]) / dspLagr;
            depthImg2 = imresize(depthImg2, [nR, nC]) / dspLagr;
            if methodID == 7 || strcmp(methodID, 'lagrSmooth'),
                maxMul = max((nA+0.5) * dA - 0.5, 0.5);
                depthImg1 = horizontalScanSmooth(depthImg1, maxMul, 3);
                depthImg2 = horizontalScanSmooth(depthImg2, -maxMul, 3);
            end
            if runParam.rho < 1e-3,
                for vId = 1:2*nA,
                    disp(sprintf('View %d: ', vId));
                    mul = abs(nA+0.5-vId) * dA - 0.5;
                    out_file = sprintf('%s/%03d_%03d.jpg', out_folder, iImg, vId);
                    if noRecomp && exist(out_file, 'file') == 2,
                        disp('not recomputing');
                        continue;
                    end
                    if mul < 0,
                        if vId > nA,
                            mul = -mul-1;
                        end
                        [outimg1, mask1] = warpDisp(img1, depthImgNonFilt1, mul, 2, true);
                        [outimg2, mask2] = warpDisp(img2, depthImgNonFilt2, mul+1, 2, true);
                        mask1 = maskFilter(mask1, 3, 10);
                        mask2 = maskFilter(mask2, 3, 10);
                        mix = (1+mul) + mask1 * (-1-mul) + mask2 * (-mul);
                        mix = repmat(mix, [1, 1, size(img1, 3)]);
                        outimg = outimg1 .* mix + outimg2 .* (1-mix);
                    elseif vId <= nA,
                        outimg = warpDisp(img1, depthImg1, mul, 2, true);
                    else,
                        outimg = warpDisp(img2, depthImg2, -mul, 2, true);
                    end
                    imwrite(imresize(outimg, [nR, nCOrig]), out_file,'jpg');
                end
            else,
                out_file = sprintf('%s/%03d_%03d.jpg', out_folder, iImg, 1);
                if noRecomp && exist(out_file, 'file') == 2,
                    disp('not recomputing');
                    continue;
                end
                weight = zeros(2*nA);
                sumImg = zeros([nR, nC, nCol, 2*nA]);
                step = 0.2;
                for vId = 1-floor(2*rho/step)*step:step:2*nA+2*rho,
                    flag = 0;
                    for uId = 1:2*nA,
                        if abs(uId - vId) < 2.1 * rho,
                            flag = 1;
                        end
                    end
                    if ~flag,
                        continue;
                    end
                    disp(sprintf('Vid %.3f: ', vId));
                    mul = abs(nA+0.5-vId) * dA - 0.5;
                    if mul < -0.05,
                        if vId > nA + 0.5,
                            mul = -mul-1;
                        end
                        [outimg1, mask1] = warpDisp(img1, depthImgNonFilt1, mul, 2, true);
                        [outimg2, mask2] = warpDisp(img2, depthImgNonFilt2, mul+1, 2, true);
                        mask1 = maskFilter(mask1, 3, 10);
                        mask2 = maskFilter(mask2, 3, 10);
                        mix = (1+mul) + mask1 * (-1-mul) + mask2 * (-mul);
                        mix = repmat(mix, [1, 1, size(img1, 3)]);
                        outimg = outimg1 .* mix + outimg2 .* (1-mix);
                    elseif vId < nA + 0.5,
                        outimg = warpDisp(img1, depthImg1, mul, 2, true);
                    else,
                        outimg = warpDisp(img2, depthImg2, -mul, 2, true);
                    end
                    for uId = 1:2*nA,
                        if abs(uId - vId) < 2.1 * rho,
                            w = exp(-0.5 * (uId - vId)^2 / (rho^2));
                            sumImg(:, :, :, uId) = sumImg(:, :, :, uId) + outimg * w;
                            weight(uId) = weight(uId) + w;
                        end
                    end
                end
                for vId = 1:2*nA,
                    out_file = sprintf('%s/%03d_%03d.jpg', out_folder, iImg, vId);
                    imwrite(imresize(sumImg(:, :, :, vId) / weight(vId), [nR, nCOrig]), out_file,'jpg');
                end
            end
        elseif methodID == 3 || strcmp(methodID, 'eul2'),
            outimg = eularian2D(img1, img2, testParam, runParam);
            nR = size(img1, 1);
            nC = size(img1, 2);

            for m=1:size(outimg, 4);
                imwrite(imresize(outimg(:, :, :, m), [nR, nCOrig]), sprintf('%s/%03d_%03d.jpg', out_folder, iImg, m),'jpg');
            end
        elseif methodID == 4 || strcmp(methodID, 'lagrWav'),
            mem = 0;
            nA = runParam.nA;
            nR = size(img1, 1);
            nC = size(img1, 2);
            for vId = 1:2*nA,
                vId
                [synth, mem, phsCmpO] = eulaOnLagrComp(reflectEx(img1), reflectEx(img2), reflectEx(depthImg1), reflectEx(depthImg2), testParam, runParam, mem, vId, true, phsCmpO);

                imwrite(imresize(synth(:, 1:nC, :), [nR, nCOrig]), sprintf('%s/%03d_%03d.jpg', out_folder, iImg, vId),'jpg');
            end
        elseif methodID == 5 || strcmp(methodID, 'eulaProp'),
            mem = 0;
            nA = runParam.nA;
            nR = size(img1, 1);
            nC = size(img1, 2);

            for vId = 1:2*nA,
                disp(sprintf('View %d: ', vId));
                out_file = sprintf('%s/%03d_%03d.jpg', out_folder, iImg, vId);
                if noRecomp && exist(out_file, 'file') == 2,
                    disp('not recomputing');
                    continue;
                end
                [synth, mem, phsCmpO] = eulaOnLagrProp(reflectEx(img1), reflectEx(img2), 0, 0, testParam, runParam, mem, vId, false, phsCmpO);
                nCRef = size(synth, 2);
                nEx = (nCRef - nC)/2;
                imwrite(imresize(synth(:, nEx+(1:nC), :), [nR, nCOrig]), out_file,'jpg');
            end
        elseif methodID == 6 || strcmp(methodID, 'eulaPropHf'),
            mem = 0;
            nA = runParam.nA;
            nR = size(img1, 1);
            nC = size(img1, 2);

            for vId = 1:2*nA,
                disp(sprintf('View %d: ', vId));
                out_file = sprintf('%s/%03d_%03d.jpg', out_folder, iImg, vId);
                if noRecomp && exist(out_file, 'file') == 2,
                    disp('not recomputing');
                    continue;
                end
                [synth, mem, phsCmpO] = eulaOnLagrProp(reflectEx(img1), reflectEx(img2), reflectEx(depthImg1), reflectEx(depthImg2), testParam, runParam, mem, vId, false, phsCmpO);
                nCRef = size(synth, 2);
                nEx = (nCRef - nC)/2;
                imwrite(imresize(synth(:, nEx+(1:nC), :), [nR, nCOrig]), out_file,'jpg');
            end
        else,
            error('unknown method');
        end

        disp(sprintf('running time: %f s', cputime() - t0));
    end
end

function filtered = maskFilter(mask, sig, hw)
    nR = size(mask, 1);
    nC = size(mask, 2);
    filtered = zeros(nR, nC);
    maskEx = zeros(nR+2*hw, nC+2*hw);
    maskEx(hw+(1:nR), hw+(1:nC)) = mask;
    for dR = -hw:hw,
        for dC = -hw:hw,
            filtered = max(filtered, maskEx(hw+dR+(1:nR), hw+dC+(1:nC)) * exp(-(dR^2+dC^2)/(2*sig^2)));
        end
    end
end