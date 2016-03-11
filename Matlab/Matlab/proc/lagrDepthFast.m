function [depthL, depthR, depthMatchL, depthMatchR, nextVals, debugV] = lagrDepthFast(Il, Ir, dStart, dStep, nD, dsR, dsC, kR, kC, kRM, kCM, useColorGuide, prevVals, alpha)
    if nargin < 13,
        prevVals = struct();
        alpha = 0;
    end
    if ~isfield(prevVals, 'L'),
        prevVals.L = cell(nD, 1);
        prevVals.R = cell(nD, 1);
        prevVals.L0 = cell(nD, 1);
        prevVals.R0 = cell(nD, 1);
    end
    nextVals.L = cell(nD, 1);
    nextVals.R = cell(nD, 1);
    nextVals.L0 = cell(nD, 1);
    nextVals.R0 = cell(nD, 1);
    debugV = struct();

    thresColor = 7/255;     % \tau_1 in eq. (5)
    thresGrad = 3/255;      % \tau_2 in eq. (5)
    gamma = 0.11;           % (1- \alpha) in eq. (5)
    threshBorder = 3/255;   % some threshold for border pixels


    [m,n,c] = size(Il);

    Il_g = rgb2gray(Il);
    Ir_g = rgb2gray(Ir);

    fx_l = gradient(Il_g);
    fx_r = gradient(Ir_g);
    fx_l = fx_l+0.5; % To get a range of values between 0 to 1
    fx_r = fx_r+0.5; % To get a range of values between 0 to 1

    ms = floor(m/dsR);
    ns = floor(n/dsC);
    Il_small = imresize(Il, [ms, ns]);
    Ir_small = imresize(Ir, [ms, ns]);
    Il_g_small = imresize(Il_g, [ms, ns]);
    Ir_g_small = imresize(Ir_g, [ms, ns]);

    depthMatchL = zeros(ms, ns);
    depthMatchR = zeros(ms, ns);
    costL = 1e4*ones(ms, ns);
    costR = 1e4*ones(ms, ns);
    costBestPrvL = zeros(ms, ns);
    costBestPrvR = zeros(ms, ns);
    costBestNxtL = zeros(ms, ns);
    costBestNxtR = zeros(ms, ns);
    prvCostL = zeros(ms, ns);
    prvCostR = zeros(ms, ns);

    for iD = 1:nD,
        dCur = dStart + dStep * (iD - 1);
        errL = gamma * min(sum(abs(extShift(Ir, -dCur) - Il), 3) / 3, thresColor) + ...
               (1 - gamma) * min(abs(extShift(fx_r, -dCur) - fx_l), thresGrad);
        errR = gamma * min(sum(abs(extShift(Il, dCur) - Ir), 3) / 3, thresColor) + ...
               (1 - gamma) * min(abs(extShift(fx_l, dCur) - fx_r), thresGrad);
        errL = (errL + 1e-3) * (1 + 0.1*(dCur / (dStep * nD)) ^ 2);
        errR = (errR + 1e-3) * (1 + 0.1*(dCur / (dStep * nD)) ^ 2);
        if useColorGuide,
            [errL, nextVals.L0{iD}] = guidedfilter_resize(errL, 1, Il, kR, kC, Il_small, Il_small, prevVals.L0{iD}, alpha);
            [errR, nextVals.R0{iD}] = guidedfilter_resize(errR, 1, Ir, kR, kC, Ir_small, Ir_small, prevVals.R0{iD}, alpha);
        else,
            [errL, nextVals.L0{iD}] = guidedfilter_resize(errL, 1, Il_g, kR, kC, Il_g_small, Il_g_small, prevVals.L0{iD}, alpha);
            [errR, nextVals.R0{iD}] = guidedfilter_resize(errR, 1, Ir_g, kR, kC, Ir_g_small, Ir_g_small, prevVals.R0{iD}, alpha);
        end
        updateL = errL < costL;
        updateR = errR < costR;
        costL(updateL) = errL(updateL);
        costR(updateR) = errR(updateR);
        depthMatchL(updateL) = iD;
        depthMatchR(updateR) = iD;

        prvUpdateL = depthMatchL == iD-1;
        prvUpdateR = depthMatchR == iD-1;
        costBestNxtL(prvUpdateL) = errL(prvUpdateL);
        costBestNxtR(prvUpdateR) = errR(prvUpdateR);

        costBestPrvL(updateL) = prvCostL(updateL);
        costBestPrvR(updateR) = prvCostR(updateR);
        prvCostL = errL;
        prvCostR = errR;
    end
    updateFracL = (depthMatchL ~= 1) & (depthMatchL ~= nD);
    updateFracR = (depthMatchR ~= 1) & (depthMatchR ~= nD);
    depthFracL = (costBestPrvL - costBestNxtL) ./ max(costBestPrvL - costL, costBestNxtL - costL);
    depthFracR = (costBestPrvR - costBestNxtR) ./ max(costBestPrvR - costR, costBestNxtR - costR);
    depthMatchL(updateFracL) = depthMatchL(updateFracL) + 0.5 * depthFracL(updateFracL);
    depthMatchR(updateFracR) = depthMatchR(updateFracR) + 0.5 * depthFracR(updateFracR);

    debugV.depthOrigL = depthMatchL;
    debugV.depthOrigR = depthMatchR;

    % Left-right consistency check
    Y = repmat((1:ms)', [1 ns]);
    X = repmat(1:ns, [ms 1]) + round(((depthMatchL-1) * dStep + dStart) / dsC);
    X(X<1) = 1;
    X(X>ns) = ns;
    indices = sub2ind([ms,ns],Y,X);
    checkLRmaskL = abs(depthMatchL - depthMatchR(indices)) < 0.5;

    Y = repmat((1:ms)', [1 ns]);
    X = repmat(1:ns, [ms 1]) - round(((depthMatchR-1) * dStep + dStart) / dsC);
    X(X<1) = 1;
    X(X>ns) = ns;
    indices = sub2ind([ms,ns],Y,X);
    checkLRmaskR = abs(depthMatchR - depthMatchL(indices)) < 0.5;

    debugV.maskL = checkLRmaskL;
    debugV.maskR = checkLRmaskR;

    depthMatchL = fillReference(depthMatchL, checkLRmaskL);
    depthMatchR = fillReference(depthMatchR, checkLRmaskR);

    sumL = zeros(m, n);
    sumR = zeros(m, n);
    depthL = zeros(m, n);
    depthR = zeros(m, n);
    for iD = 1:nD,
        valL = max(1 - abs(depthMatchL - iD), 0);
        valR = max(1 - abs(depthMatchR - iD), 0);

        minW = 0.1;
        if useColorGuide,
            [valL] = guidedfilter_resize(valL, max(checkLRmaskL, minW), Il_small, kRM, kCM, Il_small, Il);
            [valR] = guidedfilter_resize(valR, max(checkLRmaskR, minW), Ir_small, kRM, kCM, Ir_small, Ir);
        else,
            [valL] = guidedfilter_resize(valL, max(checkLRmaskL, minW), Il_g_small, kRM, kCM, Il_g_small, Il_g);
            [valR] = guidedfilter_resize(valR, max(checkLRmaskR, minW), Ir_g_small, kRM, kCM, Ir_g_small, Ir_g);
        end

        updateL = (sumL < 0.5) & (sumL + valL > 0.5);
        updateR = (sumR < 0.5) & (sumR + valR > 0.5);
        depthMedianL = iD - 0.5 + (0.5 - sumL) ./ valL;
        depthMedianR = iD - 0.5 + (0.5 - sumR) ./ valR;
        depthL(updateL) = depthMedianL(updateL);
        depthR(updateR) = depthMedianR(updateR);
        sumL = sumL + valL;
        sumR = sumR + valR;
    end
end


function res = fillReference(val, mask) 
    [m, n] = size(val);
    res = zeros(m, n);
    lastV = 1e4 * ones(m, 1);
    for iC = 1:n,
        update = mask(:, iC) ~= 0;
        valSlice = val(:, iC);
        lastV(update) = valSlice(update);
        res(:, iC) = lastV;
    end
    lastV = 1e4 * ones(m, 1);
    for iC = n:-1:1,
        update = mask(:, iC) ~= 0;
        valSlice = val(:, iC);
        lastV(update) = valSlice(update);
        res(:, iC) = min(res(:, iC), lastV);
    end
end