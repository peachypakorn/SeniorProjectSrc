tests = [];
tests = [tests; 1, 18, 1, 1, 1];


view_pos = 20;
view_dist = 400;

close all;
warning('off', 'images:initSize:adjustingMag');
for nT = 1:size(tests, 1),
    [fileLoc, testParam] = getTestCaseComp(tests(nT,1));
    runParam = getRunParamComp(tests(nT,1), tests(nT,2));
    getSynthView(fileLoc, testParam, runParam, tests(nT,3), tests(nT, 4), tests(nT, 5), view_pos, view_dist);
end

