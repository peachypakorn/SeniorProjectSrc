test = [1, 4, 1, 1, 1];

[fileLoc, testParam] = getTestCase(test(1));
runParam = getRunParam(test(1), test(2));
imgs = getOutput(fileLoc, testParam, runParam, test(3), test(4), test(5));
nR = size(imgs, 1);
nC = size(imgs, 2);
nImg = size(imgs, 5);
mul = 3;
imseq = zeros(nR*mul, nC*mul, 3, nImg+1);
mask = zeros(nR*mul, nC*mul, 3);
for iR = 1:nR*mul,
    for iC = 1:nC*mul,
        for iCol = 1:3,
            mask(iR, iC, iCol) = mod(iR-1 + (iC-1) * 3 + iCol-1, 8) + 1;
        end
    end
end
for iImg = 1:nImg,
    imtmp = zeros(nR*mul, nC*mul, 3);
    for iA = 1:8,
        it = imresize(imgs(:, :, :, iA, iImg), mul);;
        imtmp(mask == iA) = it(mask == iA);
    end
    imseq(:, :, :, iImg) = imtmp;
end
imseq(:, :, :, nImg+1) = imseq(:, :, :, 1);
player = implay(uint8(imseq));