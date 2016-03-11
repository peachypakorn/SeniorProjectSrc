function [filteredVal, nextVals] = guidedfilter_resize(val, weight, guide_input, kR, kC, guide_box, guide_output, prevVals, alpha)
    if nargin < 8
        prevVals = struct();
        alpha = 0;
    end

    if size(guide_input, 3) > 3,
        disp('Error: #channels > 3');
        filteredVal = (guidedfilter_resize(val, guide_input(:, :, 1), kR, kC, guide_box(:, :, 1), guide_output(:, :, 1)) + ...
                       guidedfilter_resize(val, guide_input(:, :, 2), kR, kC, guide_box(:, :, 2), guide_output(:, :, 2)) + ...
                       guidedfilter_resize(val, guide_input(:, :, 3), kR, kC, guide_box(:, :, 3), guide_output(:, :, 3))) / 3;
        return
    end
    [ri, ci, chi] = size(guide_input);
    [r, c, ch] = size(guide_box);
    [ro, co, cho] = size(guide_output);
    if numel(weight) == 1,
        weight = ones(ri, ci);
    end
    if ~isfield(prevVals, 'covVal'),
        prevVals.covVal = zeros([r, c, ch+1]);
        prevVals.covGuide = zeros([r, c, ch, ch]);
        prevVals.meanGuide = zeros([r, c, ch]);
        prevVals.weightGuide = zeros([r, c]);
        alpha = 0;
    end
    nextVals = struct();

    covVal = zeros([r, c, ch+1]);
    for iCh = 1:ch,
        covVal(:, :, iCh) = imresize(weight .* val .* guide_input(:, :, iCh), [r, c], 'bilinear');
    end
    covVal(:, :, ch+1) = imresize(weight .* val, [r, c], 'bilinear');
    covVal = boxfilter(covVal, [kR, kC]);
    meanGuide = boxfilter(imresize(repmat(weight, [1, 1, ch]) .* guide_input, [r, c], 'bilinear'), [kR, kC]);

    covGuide = zeros([r, c, ch, ch]);
    for iCh = 1:ch,
        for jCh = iCh:ch,
            covGuide(:, :, iCh, jCh) = boxfilter(imresize(weight .* guide_input(:, :, iCh) ...
                    .* guide_input(:, :, jCh), [r, c], 'bilinear'), [kR, kC]);
        end
    end
    weightGuide = boxfilter(imresize(weight, [r, c], 'bilinear'), [kR, kC]);

    covVal = alpha * prevVals.covVal + (1-alpha) * covVal;
    covGuide = alpha * prevVals.covGuide + (1-alpha) * covGuide;
    meanGuide = alpha * prevVals.meanGuide + (1-alpha) * meanGuide;
    weightGuide = alpha * prevVals.weightGuide + (1-alpha) * weightGuide;
    if nargin >= 8
        nextVals.covVal = covVal;
        nextVals.covGuide = covGuide;
        nextVals.meanGuide = meanGuide;
        nextVals.weightGuide = weightGuide;
    end
    covGuide = covGuide ./ repmat(weightGuide, [1, 1, ch, ch]);
    meanGuide = meanGuide ./ repmat(weightGuide, [1, 1, ch]);
    covVal = covVal ./ repmat(weightGuide, [1, 1, ch+1]);

    covVal(:, :, 1:ch) = covVal(:, :, 1:ch) - meanGuide .* repmat(covVal(:, :, ch+1), [1, 1, ch]);

    for iCh = 1:ch,
        for jCh = iCh:ch,
            covGuide(:, :, iCh, jCh) = covGuide(:, :, iCh, jCh) - meanGuide(:, :, iCh) .* meanGuide(:, :, jCh);
        end
        covGuide(:, :, iCh, iCh) = covGuide(:, :, iCh, iCh) + 1e-3;
    end
    for iCh = 1:ch,
        if 0,
            % Not needed since we only use terms with iCh <= jCh
            for jCh = 1:iCh-1,
                covGuide(:, :, iCh, jCh) = covGuide(:, :, jCh, iCh);
            end
        end
    end

    covInv = zeros([r, c, ch, ch]);
    covMul = zeros([r, c, ch]);
    if ch == 1,
        covDet = covGuide(:, :, 1, 1);
        covDet = max(covDet, 1e-3);
        covInv(:, :, 1, 1) = 1 ./ covDet;
        covMul(:, :, 1) = covVal(:, :, 1) .* covInv(:, :, 1, 1);
    elseif ch == 2,
        covDet = covGuide(:, :, 1, 1) .* covGuide(:, :, 2, 2) - covGuide(:, :, 1, 2) .^ 2;
        covDet = max(covDet, 1e-5);
        covInv(:, :, 1, 1) = covGuide(:, :, 2, 2) ./ covDet;
        covInv(:, :, 1, 2) = -covGuide(:, :, 1, 2) ./ covDet;
        covInv(:, :, 2, 2) = covGuide(:, :, 1, 1) ./ covDet;
        covMul(:, :, 1) = covVal(:, :, 1) .* covInv(:, :, 1, 1) + covVal(:, :, 2) .* covInv(:, :, 1, 2);
        covMul(:, :, 2) = covVal(:, :, 1) .* covInv(:, :, 1, 2) + covVal(:, :, 2) .* covInv(:, :, 2, 2);
    else,
        covInv(:, :, 1, 1) = covGuide(:, :, 2, 2) .* covGuide(:, :, 3, 3) - covGuide(:, :, 2, 3) .^ 2;
        covInv(:, :, 1, 2) = covGuide(:, :, 2, 3) .* covGuide(:, :, 1, 3) - covGuide(:, :, 1, 2) .* covGuide(:, :, 3, 3);
        covInv(:, :, 1, 3) = covGuide(:, :, 1, 2) .* covGuide(:, :, 2, 3) - covGuide(:, :, 1, 3) .* covGuide(:, :, 2, 2);
        covInv(:, :, 2, 2) = covGuide(:, :, 1, 1) .* covGuide(:, :, 3, 3) - covGuide(:, :, 1, 3) .^ 2;
        covInv(:, :, 2, 3) = covGuide(:, :, 1, 2) .* covGuide(:, :, 1, 3) - covGuide(:, :, 2, 3) .* covGuide(:, :, 1, 1);
        covInv(:, :, 3, 3) = covGuide(:, :, 1, 1) .* covGuide(:, :, 2, 2) - covGuide(:, :, 1, 2) .^ 2;
        covDet = sum(covGuide(:, :, 1, :) .* covInv(:, :, 1, :), 4);
        %covDet = max(covDet, 1e-7);
        covMul(:, :, 1) = covVal(:, :, 1) .* covInv(:, :, 1, 1) + covVal(:, :, 2) .* covInv(:, :, 1, 2) + covVal(:, :, 3) .* covInv(:, :, 1, 3);
        covMul(:, :, 2) = covVal(:, :, 1) .* covInv(:, :, 1, 2) + covVal(:, :, 2) .* covInv(:, :, 2, 2) + covVal(:, :, 3) .* covInv(:, :, 2, 3);
        covMul(:, :, 3) = covVal(:, :, 1) .* covInv(:, :, 1, 3) + covVal(:, :, 2) .* covInv(:, :, 2, 3) + covVal(:, :, 3) .* covInv(:, :, 3, 3);
        covMul = covMul ./ repmat(covDet, [1, 1, 3]);
    end

    covVal(:, :, 1:ch) = covMul;
    covVal(:, :, ch+1) = boxfilter(imresize(val, [r, c], 'bilinear'), [kR, kC]) - ...
            sum(covVal(:, :, 1:ch) .* boxfilter(guide_box(:, :, :), [kR, kC]), 3);
    covVal = boxfilter(covVal, [kR, kC]);

    filteredVal = sum(imresize(covVal(:, :, 1:ch), [ro, co], 'bilinear') .* guide_output(:, :, :), 3) + imresize(covVal(:, :, ch+1), [ro, co], 'bilinear');
end