function getSynthView(fileID, testParam, runParam, methodID, startInd, dsId, view_pos, view_dist)
    in_folder = sprintf('result/compare/%s/%s/%s-%s', fileID.testName, ...
            runParam.runName, getMethodNameComp(methodID), getDepthSelect(dsId));
    out_folder = sprintf('%s/synth', in_folder);
    mkdir(out_folder);

    for iImg = startInd:fileID.nImg,
        disp(sprintf('%s-%s-%s-%d-%s', fileID.testName, runParam.runName, getMethodNameComp(methodID), iImg, getDepthSelect(dsId)));
        t1 = cputime();
        for iA = 1:2*runParam.nA;
            filename = sprintf('%s/%03d_%03d.jpg', in_folder, iImg, iA);
            img = imread(filename);
            if iA == 1,
                imgs1 = zeros([size(img), 2*runParam.nA]);
            end
            imgs1(:, :, :, iA) = double(img)/255;
        end

        [outimgL, outimgR] = viewSimulate(imgs1, view_pos, view_dist);
        imwrite(outimgL, sprintf('%s/v_%03d_%03d_L.jpg', out_folder, view_dist, view_pos), 'jpg');
        imwrite(outimgR, sprintf('%s/v_%03d_%03d_R.jpg', out_folder, view_dist, view_pos), 'jpg');
        outimgBWL = outimgL(:, :, 1) * 0.3 + outimgL(:, :, 2) * 0.6 + outimgL(:, :, 3) * 0.1;
        outimgBWR = outimgR(:, :, 1) * 0.3 + outimgR(:, :, 2) * 0.6 + outimgR(:, :, 3) * 0.1;
        [depth1, depth2] = getCCDepth(outimgBWL, outimgBWR, testParam.maxDepth);
        imwrite(depth1 / (2*testParam.maxDepth) + 0.5, ...
                sprintf('%s/d_%03d_%03d_L.jpg', out_folder, view_dist, view_pos), 'jpg');
        imwrite(depth2 / (2*testParam.maxDepth) + 0.5, ...
                sprintf('%s/d_%03d_%03d_R.jpg', out_folder, view_dist, view_pos), 'jpg');
        disp(sprintf('time: %.3f s', cputime() - t1));
    end
end

