tests = [25, 12, 5, 5];

for nT = 1:size(tests, 1),
    [fileLoc, testParam] = getTestCase(tests(nT,1));
    runParam = getRunParam(tests(nT,1), tests(nT,2));
    runTest(fileLoc, testParam, runParam, tests(nT,3), tests(nT, 4));
end
