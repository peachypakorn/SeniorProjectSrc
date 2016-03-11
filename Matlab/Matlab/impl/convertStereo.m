tests = [7, 4, 1, 1, 12; 6, 4, 5, 1, 12];

for tId = 1:size(tests, 1),
    [fileLoc, testParam] = getTestCase(tests(tId, 1));
    runParam = getRunParam(tests(tId, 1), tests(tId, 2));
    imgs = getOutput(fileLoc, testParam, runParam, tests(tId, 3), tests(tId, 4), tests(tId, 5));
    nR = size(imgs, 1);
    nC = size(imgs, 2);
    nImg = size(imgs, 5);
    mul = min(3840 / nC, 2160 / nR);
    mask = permute(zeros(size(imresize(imgs(:, :, :, 1, 1), mul))), [2, 1, 3]);
    [CC, RR, VV] = meshgrid(1:size(mask, 1), 1:size(mask, 2), 1:3);
    mask = mod(RR-1 + (CC-1) * 3 + (3-VV), 8) + 1;
    out_folder = sprintf('result/viewing/%s/%s/%s/stereo', fileLoc.testName, runParam.runName, getMethodName(tests(tId, 4)));
    mkdir(out_folder)
    for iImg = 1:nImg,
        disp(sprintf('%s-%s-%s-%d', fileLoc.testName, runParam.runName, getMethodName(tests(tId, 3)), iImg));
        imtmp = zeros(size(mask));
        for iA = 1:8,
            it = imresize(imgs(:, :, :, iA, iImg), mul);;
            imtmp(mask == iA) = it(mask == iA);
        end
	    imwrite(uint8(imtmp), sprintf('%s/%03d_stereo.png', out_folder, iImg + tests(tId, 4) - 1), 'png');
    end
    clear mask;
    clear imtmp;
    clear CC;
    clear RR;
    clear VV;
    clear imgs;
end
