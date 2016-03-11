%cases = (7:14)';
%cases = (21:24)';
%cases = 1;
%cases = [1, 41:43, 18:21, 11:13]';
tests = [];

cases = [1];
caseN = length(cases);
tests = [tests; cases, 25*ones(caseN, 1), ones(caseN, 1), ones(caseN, 1), ones(caseN, 1)];
%tests = [tests; cases, ones(caseN, 1), 2*ones(caseN, 1), ones(caseN, 1), 5*ones(caseN, 1)];
%tests = [tests; cases, ones(caseN, 1), 2*ones(caseN, 1), ones(caseN, 1), 2*ones(caseN, 1)];
%tests = [tests; cases, ones(caseN, 1), ones(caseN, 1), ones(caseN, 1), 4*ones(caseN, 1)];
%tests = [tests; cases, ones(caseN, 1), ones(caseN, 1), ones(caseN, 1), ones(caseN, 1)];
%tests = [tests; 51, 17, 1, 1, 1];

noRecomp = false;

if true,
    depthRecompParam = struct();
    depthRecompParam.recompLowD = false;
    depthRecompParam.recompMedi = false;
    depthRecompParam.recompUpD = false;
    depthRecompParam.recompImpD = false;
    depthRecompParam.recompGtD = false;
    depthRecompParam.recompHafD = false;
    depthRecompParam.nocompHafD = true;
    depthRecompParam.compLowImpD = false;
else,
    depthRecompParam = 0;
end
close all;
warning('off', 'images:initSize:adjustingMag');
for nT = 1:size(tests, 1),
    [fileLoc, testParam] = getTestCaseComp(tests(nT,1));
    runParam = getRunParamComp(tests(nT,1), tests(nT,2));
    runTestComp(fileLoc, testParam, runParam, tests(nT,3), tests(nT, 4), tests(nT, 5), depthRecompParam, noRecomp);
end
