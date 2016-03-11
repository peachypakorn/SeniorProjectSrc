tests = [5, 11, 5, 1];

for nT = 1:size(tests, 1),
    [fileLoc, testParam] = getTestCase(tests(nT,1));
    runParam = getRunParam(tests(nT,1), tests(nT,2));
    runTest(fileLoc, testParam, runParam, tests(nT,3), tests(nT, 4));
end
