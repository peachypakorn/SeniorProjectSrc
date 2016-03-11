nC = 256;
sp1 = 100;
l = 20;
d = 12;

img1 = [zeros(10, sp1, 1), ones(10, l, 1), zeros(10, nC-sp1-l, 1)];
img2 = [zeros(10, sp1+d, 1), ones(10, l, 1), zeros(10, nC-sp1-l-d, 1)];
depth1 = [zeros(10, sp1, 1), -d*ones(10, l, 1), zeros(10, nC-sp1-l, 1)];
depth2 = [zeros(10, sp1+d, 1), -d*ones(10, l, 1), zeros(10, nC-sp1-l-d, 1)];
%img1 = exp(-4*((1:nC) - sp1).^2 / l^2);
%img2 = exp(-4*((1:nC) - sp1-d).^2 / l^2);
%depth1 = [zeros(10, sp1-l, 1), -d*ones(10, 2*l, 1), zeros(10, nC-sp1-l, 1)];
%depth2 = [zeros(10, sp1-l+d, 1), -d*ones(10, 2*l, 1), zeros(10, nC-sp1-l-d, 1)];
%depth1 = zeros(nC, 1);
%depth2 = zeros(nC, 1);

testParam = struct();
testParam.maxDepth = 0;
runParam = struct();
runParam.rho = 0.0;
runParam.ht = 6;
runParam.nA = 1;
runParam.dA = 7;
runParam.eurLagrResize = 0.25;

[out, vis] = eulaOnLagrVis(img1, img2, depth1, depth2, testParam, runParam, false);

figure;
view = 1;
set(gca,'Visible','off');
for iLay = runParam.ht+2:-1:1,
    h = subplot(runParam.ht+2, 2, (runParam.ht+3 - iLay)*2-1);
    p = get(h, 'pos');
    p(4) = 0.1;
    set(h, 'pos', p);
    hold on;
    nL = length(vis.pyr1{iLay});
    nL2 = length(vis.mov{1, 1}.sig{iLay});
    xlim([0, nL]);
    xRange = 0.5 + (0:nL-1);
    xRange2 = 0.5 + (nL / nL2) * (0:nL2-1);
    xRangeC = (0.5 + (0:nC-1)) * nL / nC;
    plot(xRange, real(vis.pyr1{iLay}), '*', 'MarkerSize', 2);
    plot(xRangeC, vis.layer1{iLay});
    %plot(xRangeC, vis.accum1{iLay});
    %stem(vis.mov{view, 1}.mov{iLay} * nL / nL2, real(vis.pyr1{iLay}), 'k.', 'MarkerSize', 2);
    %plot(xRange2, real(vis.mov{view, 1}.sig{iLay}), 'k*', 'MarkerSize', 2);
    plot(xRange, real(vis.mov{view, 1}.sigFilt{iLay}), 'r*', 'MarkerSize', 2);
    plot(xRangeC, vis.mov{view, 1}.layer{iLay}, 'red');
    %plot(xRangeC, vis.mov{view, 1}.recon{iLay}, 'red');
end

for iLay = runParam.ht+2:-1:1,
    h = subplot(runParam.ht+2, 2, (runParam.ht+3 - iLay)*2);
    p = get(h, 'pos');
    p(4) = 0.1;
    set(h, 'pos', p);
    hold on;
    nL = length(vis.pyr1{iLay});
    nL2 = length(vis.mov{1, 1}.sig{iLay});
    xlim([0, nL]);
    xRangeC = (0.5 + (0:nC-1)) * nL / nC;
    plot(xRangeC, vis.accum1{iLay});
    plot(xRangeC, vis.mov{view, 1}.recon{iLay}, 'red');
end
