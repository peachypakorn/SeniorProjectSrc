function runTest(fileID, testParam, runParam, methodID, startInd)
    out_folder = sprintf('result/all/%s/%s/%s', fileID.testName, runParam.runName, getMethodName(methodID));
    mkdir(out_folder);

    phsCmpO = 0;
    for iImg = startInd:fileID.nImg,
        disp(sprintf('%s-%s-%s-%d', fileID.testName, runParam.runName, getMethodName(methodID), iImg));
        t0 = cputime();
        if strcmp(fileID.inpType, 'separate'),
            img1 = double(imread(sprintf('%s%s', fileID.inpPrefix, fileID.inpImg1{iImg}))) / 255;
            img2 = double(imread(sprintf('%s%s', fileID.inpPrefix, fileID.inpImg2{iImg}))) / 255;
        elseif strcmp(fileID.inpType, 'combined'),
            img = double(imread(sprintf('%s%s', fileID.inpPrefix, fileID.inpImg{iImg}))) / 255;
            nC = size(img, 2);
            img1 = img(:, 1:nC/2, :);
            img2 = img(:, nC/2+1:nC, :);
        elseif strcmp(fileID.inpType, 'combinedInv'),
            img = double(imread(sprintf('%s%s', fileID.inpPrefix, fileID.inpImg{iImg}))) / 255;
            nC = size(img, 2);
            img2 = img(:, 1:nC/2, :);
            img1 = img(:, nC/2+1:nC, :);
        else,
            error('unknown inpType');
        end
        %img1 = img1 .^ 2.2;
        %img2 = img2 .^ 2.2;

        if isfield(fileID, 'resize'),
            img1 = imresize(img1, fileID.resize);
            img2 = imresize(img2, fileID.resize);
        end
        depthImg1 = 0;
        depthImg2 = 0;
        if isfield(fileID, 'depthImg1'),
            if exist(fileID.depthImg1{iImg}, 'file') == 2,
                depthImg1 = round((double(imread(fileID.depthImg1{iImg})) / 255 - 0.5) * 2 * fileID.maxDepth);
            end
            if exist(fileID.depthImg2{iImg}, 'file') == 2,
                depthImg2 = round((double(imread(fileID.depthImg2{iImg})) / 255 - 0.5) * 2 * fileID.maxDepth);
            end
        end
        if (isscalar(depthImg1) || isscalar(depthImg2)) && methodID ~= 3,
            disp('no depth image found');
        end

        if methodID == 1 || strcmp(methodID, 'eulaOnLagr'),
            [outimg, depthL, depthR, phsCmpO] = eulaOnLagr(img1, img2, depthImg1, depthImg2, testParam, runParam, false, phsCmpO);
            if isscalar(depthImg1) || isscalar(depthImg2),
                imwrite(0.5*depthL/testParam.maxDepth + 0.5, sprintf('%s/%03d_depthL.png', out_folder, iImg),'png');
                imwrite(0.5*depthR/testParam.maxDepth + 0.5, sprintf('%s/%03d_depthR.png', out_folder, iImg),'png');
            end
            if ~isscalar(phsCmpO),
                imwrite(imresize(0.5 + phsCmpO'*runParam.nA/(4*pi), [size(img1, 1), size(img1, 2)]), sprintf('%s/%03d_phsCmp.jpg', out_folder, iImg),'jpg');
            end
        elseif methodID == 2 || strcmp(methodID, 'lagrangian'),
            [outimg, depthL, depthR] = lagrangian(img1, img2, depthImg1, depthImg2, testParam, runParam);

            imwrite(0.5*depthL/testParam.maxDepth + 0.5, sprintf('%s/%03d_depthL.png', out_folder, iImg),'png');
            imwrite(0.5*depthR/testParam.maxDepth + 0.5, sprintf('%s/%03d_depthR.png', out_folder, iImg),'png');
        elseif methodID == 3 || strcmp(methodID, 'eularian2D'),
            outimg = eularian2D(img1, img2, testParam, runParam);
        elseif methodID == 4 || strcmp(methodID, 'eulaSynth'),
            [outimg, depthL, depthR, phsCmpO] = eulaOnLagr(img1, img2, depthImg1, depthImg2, testParam, runParam, true, phsCmpO);
        elseif methodID == 5 || strcmp(methodID, 'eulaOnLagrRefl'),
            [outimg, depthL, depthR, phsCmpO] = eulaOnLagr(reflectEx(img1), reflectEx(img2), reflectEx(depthImg1), reflectEx(depthImg2), testParam, runParam, false, phsCmpO);
            nC = size(img1, 2);
            outimg = outimg(:, 1:nC, :, :);
        else,
            error('unknown method');
        end

        if runParam.runID == 6 ,
            imwrite(repmat(permute(squeeze(outimg(1, :, :, :)), [3, 1, 2]), [3, 1, 1]), sprintf('%s/%03d_lf.png', out_folder, iImg),'png');
        else
            outimg(outimg < 0) = 0;
            %outimg = outimg .^ (1/2.2);
            for m=1:size(outimg, 4);
                imwrite(outimg(:, :, :, m), sprintf('%s/%03d_%03d.jpg', out_folder, iImg, m),'jpg');
            end
        end
        disp(sprintf('running time: %f s', cputime() - t0));
    end
end


function reflect = reflectEx(img)
    if isscalar(img)
        reflect = img;
    else
        nR = size(img, 1);
        nC = size(img, 2);
        nCh = size(img, 3);
        reflect = zeros(nR, 2*nC, nCh);
        reflect(:, 1:nC, :) = img;
        reflect(:, 2*nC:-1:nC+1, :) = img;
    end
end