%cases = (7:14)';
%cases = (21:24)';
tests = [];

cases = [77]';
caseN = length(cases);
tests = [tests; cases, 29*ones(caseN, 1), ones(caseN, 1), 2*ones(caseN, 1), 9*ones(caseN, 1)];
cases = [70]';
caseN = length(cases);
tests = [tests; cases, 27*ones(caseN, 1), 3*ones(caseN, 1), 1*ones(caseN, 1), 9*ones(caseN, 1)];
tests = [tests; cases, 27*ones(caseN, 1), 7*ones(caseN, 1), 1*ones(caseN, 1), 9*ones(caseN, 1)];
tests = [tests; cases, 27*ones(caseN, 1), ones(caseN, 1), 1*ones(caseN, 1), 9*ones(caseN, 1)];

if 0
    cases = [46, 75, 45, 76, 51, 77, 50, 74, 72, 73, 64, 65, 78, 67, 68, 69, 70, 71]';
    runPs = [28, 28, 28, 27, 26, 29, 27, 27, 26, 29, 28, 28, 28, 28, 28, 28, 27, 28]';
    caseN = length(cases);
    tests = [tests; cases, runPs, 3*ones(caseN, 1), ones(caseN, 1), 9*ones(caseN, 1)];
end
%tests = [tests; cases, 28*ones(caseN, 1), ones(caseN, 1), ones(caseN, 1), 9*ones(caseN, 1)];
%tests = [tests; cases, ones(caseN, 1),ones(caseN, 1), ones(caseN, 1), 4*ones(caseN, 1)];
noRecomp = false;
if false,
    depthRecompParam = struct();
    depthRecompParam.recompLowD = false;
    depthRecompParam.recompMedi = false;
    depthRecompParam.recompUpD = false;
    depthRecompParam.recompImpD = false;
    depthRecompParam.recompGtD = true;
    depthRecompParam.recompHafD = false;
    depthRecompParam.nocompHafD = true;
    depthRecompParam.compLowImpD = false;
else,
    depthRecompParam = 0;
end
close all;
for nT = 1:size(tests, 1),
    [fileLoc, testParam] = getTestCaseComp(tests(nT,1));
    runParam = getRunParamComp(tests(nT,1), tests(nT,2));
    runTestComp(fileLoc, testParam, runParam, tests(nT,3), tests(nT, 4), tests(nT, 5), depthRecompParam, noRecomp);
end
