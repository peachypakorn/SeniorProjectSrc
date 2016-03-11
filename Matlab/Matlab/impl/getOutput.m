function outImg = getOutput(fileID, testParam, runParam, methodID, startInd, endInd)
    out_folder = sprintf('result/all/%s/%s/%s', fileID.testName, runParam.runName, getMethodName(methodID));

    for iImg = startInd:endInd,
        nA = runParam.nA * 2;
        for m=1:nA;
            im = imread(sprintf('%s/%03d_%03d.jpg', out_folder, iImg, m));
            if m == 1 && iImg == startInd,
                outImg = zeros([size(im), nA, endInd - startInd + 1]);
            end
            outImg(:, :, :, m, iImg - startInd + 1) = im;
        end
    end
end
