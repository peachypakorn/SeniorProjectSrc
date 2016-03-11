function [fileLoc, testParam] = getTestCase(testID)
    fileLoc = struct();
    testParam = struct();
    if testID == 1 || strcmp(testID, 'bbb-basic'),
        fileLoc.testName = 'bbb-basic';
        fileLoc.testID = 1;
        fileLoc.nImg = 1;
        fileLoc.inpPrefix = 'data/';
        fileLoc.inpType = 'separate';
        fileLoc.inpImg1 = {'frame1.jpg'};
        fileLoc.inpImg2 = {'frame2.jpg'};
        fileLoc.resize = 2/3;
        fileLoc.depthImg1 = {'result/all/bbb-basic/plain/lagrangian/001_depthL.png'};
        fileLoc.depthImg2 = {'result/all/bbb-basic/plain/lagrangian/001_depthR.png'};
        fileLoc.maxDepth = 12;

        testParam.maxDepth = 12;
    elseif testID == 2 || strcmp(testID, 'dracular-seq-215'),
        fileLoc.testName = 'dracular-seq-215';
        fileLoc.testID = 2;
        fileLoc.nImg = 12;
        fileLoc.inpPrefix = 'data/dracular-seq-215/';
        fileLoc.inpType = 'separate';
        fileLoc.inpImg1 = getNameArray('ud-%04d-0.jpg', 215:2:238);
        fileLoc.inpImg2 = getNameArray('ud-%04d-1.jpg', 215:2:238);
        fileLoc.resize = 2/3;
        fileLoc.depthImg1 = getNameArray('result/all/dracular-seq-215/plain/lagrangian/%03d_depthL.png', 1:12);
        fileLoc.depthImg2 = getNameArray('result/all/dracular-seq-215/plain/lagrangian/%03d_depthR.png', 1:12);
        fileLoc.maxDepth = 12;

        testParam.maxDepth = 12;
    elseif testID == 3 || strcmp(testID, 'dracular-seq-769'),
        fileLoc.testName = 'dracular-seq-769';
        fileLoc.testID = 3;
        fileLoc.nImg = 12;
        fileLoc.inpPrefix = 'data/dracular-seq-769/';
        fileLoc.inpType = 'separate';
        fileLoc.inpImg1 = getNameArray('ud-%04d-0.jpg', 769:2:792);
        fileLoc.inpImg2 = getNameArray('ud-%04d-1.jpg', 769:2:792);
        fileLoc.resize = 2/3;
        fileLoc.depthImg1 = getNameArray('result/all/dracular-seq-769/plain/lagrangian/%03d_depthL.png', 1:12);
        fileLoc.depthImg2 = getNameArray('result/all/dracular-seq-769/plain/lagrangian/%03d_depthR.png', 1:12);
        fileLoc.maxDepth = 12;

        testParam.maxDepth = 12;
    elseif testID == 4 || strcmp(testID, 'dracular-seq-969'),
        fileLoc.testName = 'dracular-seq-969';
        fileLoc.testID = 4;
        fileLoc.nImg = 12;
        fileLoc.inpPrefix = 'data/dracular-seq-969/';
        fileLoc.inpType = 'separate';
        fileLoc.inpImg1 = getNameArray('ud-%04d-0.jpg', 969:2:992);
        fileLoc.inpImg2 = getNameArray('ud-%04d-1.jpg', 969:2:992);
        fileLoc.resize = 2/3;
        fileLoc.depthImg1 = getNameArray('result/all/dracular-seq-969/plain/lagrangian/%03d_depthL.png', 1:12);
        fileLoc.depthImg2 = getNameArray('result/all/dracular-seq-969/plain/lagrangian/%03d_depthR.png', 1:12);
        fileLoc.maxDepth = 40;

        testParam.maxDepth = 40;
    elseif testID == 5 || strcmp(testID, 'airport-seq-110'),
        fileLoc.testName = 'airport-seq-110';
        fileLoc.testID = 5;
        fileLoc.nImg = 12;
        fileLoc.inpPrefix = 'data/airport-seq-110/airport-110-';
        fileLoc.inpType = 'combined';
        fileLoc.inpImg = getNameArray('ud-pa_%03d.jpg', 175:2:198);
        fileLoc.resize = 3/4;
        fileLoc.depthImg1 = getNameArray('result/all/airport-seq-110/plain/lagrangian/%03d_depthL.png', 1:12);
        fileLoc.depthImg2 = getNameArray('result/all/airport-seq-110/plain/lagrangian/%03d_depthR.png', 1:12);
        fileLoc.maxDepth = 20;

        testParam.maxDepth = 20;
    elseif testID == 6 || strcmp(testID, 'airport-seq-243'),
        fileLoc.testName = 'airport-seq-243';
        fileLoc.testID = 6;
        fileLoc.nImg = 12;
        fileLoc.inpPrefix = 'data/airport-seq-243/airport-243-';
        fileLoc.inpType = 'combined';
        fileLoc.inpImg = getNameArray('ud-pa_%03d.jpg', 21:2:44);
        fileLoc.resize = 3/4;
        fileLoc.depthImg1 = getNameArray('result/all/airport-seq-243/plain/lagrangian/%03d_depthL.png', 1:12);
        fileLoc.depthImg2 = getNameArray('result/all/airport-seq-243/plain/lagrangian/%03d_depthR.png', 1:12);
        fileLoc.maxDepth = 20;

        testParam.maxDepth = 20;
    elseif testID == 7 || strcmp(testID, 'skydiving-seq-40'),
        fileLoc.testName = 'skydiving-seq-40';
        fileLoc.testID = 7;
        fileLoc.nImg = 12;
        fileLoc.inpPrefix = 'data/skydiving-seq-40/skydiving-40-';
        fileLoc.inpType = 'combined';
        fileLoc.inpImg = getNameArray('ud-pa_%03d.jpg', 31:2:54);
        fileLoc.resize = 3/4;
        fileLoc.depthImg1 = getNameArray('result/all/skydiving-seq-40/plain/lagrangian/%03d_depthL.png', 1:12);
        fileLoc.depthImg2 = getNameArray('result/all/skydiving-seq-40/plain/lagrangian/%03d_depthR.png', 1:12);
        fileLoc.maxDepth = 20;

        testParam.maxDepth = 20;
    elseif testID == 8 || strcmp(testID, 'skydiving-seq-152'),
        fileLoc.testName = 'skydiving-seq-152';
        fileLoc.testID = 8;
        fileLoc.nImg = 12;
        fileLoc.inpPrefix = 'data/skydiving-seq-152/skydiving-152-';
        fileLoc.inpType = 'combined';
        fileLoc.inpImg = getNameArray('ud-pa_%03d.jpg', 25:2:48);
        fileLoc.resize = 3/4;
        fileLoc.depthImg1 = getNameArray('result/all/skydiving-seq-152/plain/lagrangian/%03d_depthL.png', 1:12);
        fileLoc.depthImg2 = getNameArray('result/all/skydiving-seq-152/plain/lagrangian/%03d_depthR.png', 1:12);
        fileLoc.maxDepth = 30;

        testParam.maxDepth = 30;
    elseif testID == 9 || strcmp(testID, 'sintel-seq-290'),
        fileLoc.testName = 'sintel-seq-290';
        fileLoc.testID = 9;
        fileLoc.nImg = 12;
        fileLoc.inpPrefix = 'data/sintel-seq-290/sintel-290-';
        fileLoc.inpType = 'combined';
        fileLoc.inpImg = getNameArray('ud-pa_%03d.jpg', 25:2:47);
        fileLoc.resize = 3/4;
        fileLoc.depthImg1 = getNameArray('result/all/sintel-seq-290/plain/lagrangian/%03d_depthL.png', 1:12);
        fileLoc.depthImg2 = getNameArray('result/all/sintel-seq-290/plain/lagrangian/%03d_depthR.png', 1:12);
        fileLoc.maxDepth = 15;

        testParam.maxDepth = 15;
    elseif testID == 10 || strcmp(testID, 'sintel-seq-702'),
        fileLoc.testName = 'sintel-seq-702';
        fileLoc.testID = 10;
        fileLoc.nImg = 12;
        fileLoc.inpPrefix = 'data/sintel-seq-702/sintel-702-';
        fileLoc.inpType = 'combined';
        fileLoc.inpImg = getNameArray('ud-pa_%03d.jpg', 25:2:47);
        fileLoc.resize = 3/4;
        fileLoc.depthImg1 = getNameArray('result/all/sintel-seq-702/plain/lagrangian/%03d_depthL.png', 1:12);
        fileLoc.depthImg2 = getNameArray('result/all/sintel-seq-702/plain/lagrangian/%03d_depthR.png', 1:12);
        fileLoc.maxDepth = 10;

        testParam.maxDepth = 10;
    elseif testID == 11 || strcmp(testID, 'heidelberg-seq-338'),
        fileLoc.testName = 'heidelberg-seq-338';
        fileLoc.testID = 11;
        fileLoc.nImg = 12;
        fileLoc.inpPrefix = 'data/heidelberg-seq-338/heidelberg-338-';
        fileLoc.inpType = 'combined';
        fileLoc.inpImg = getNameArray('ud_%03d.jpg', 1:2:24);
        fileLoc.resize = 3/4;
        fileLoc.depthImg1 = getNameArray('result/all/heidelberg-seq-338/plain/lagrangian/%03d_depthL.png', 1:12);
        fileLoc.depthImg2 = getNameArray('result/all/heidelberg-seq-338/plain/lagrangian/%03d_depthR.png', 1:12);
        fileLoc.maxDepth = 10;

        testParam.maxDepth = 10;
    elseif testID == 12 || strcmp(testID, 'demo-seq-579'),
        fileLoc.testName = 'demo-seq-579';
        fileLoc.testID = 12;
        fileLoc.nImg = 12;
        fileLoc.inpPrefix = 'data/demo-seq-579/';
        fileLoc.inpType = 'separate';
        fileLoc.inpImg1 = getNameArray('%04d-0.jpg', 579:2:602);
        fileLoc.inpImg2 = getNameArray('%04d-1.jpg', 579:2:602);
        fileLoc.resize = 3/4;
        fileLoc.depthImg1 = getNameArray('result/all/demo-seq-579/plain/lagrangian/%03d_depthL.png', 1:12);
        fileLoc.depthImg2 = getNameArray('result/all/demo-seq-579/plain/lagrangian/%03d_depthR.png', 1:12);
        fileLoc.maxDepth = 10;

        testParam.maxDepth = 10;
    elseif testID == 13 || strcmp(testID, 'demo-seq-721'),
        fileLoc.testName = 'demo-seq-721';
        fileLoc.testID = 13;
        fileLoc.nImg = 12;
        fileLoc.inpPrefix = 'data/demo-seq-721/';
        fileLoc.inpType = 'separate';
        fileLoc.inpImg1 = getNameArray('pa-%04d-0.jpg', 721:2:744);
        fileLoc.inpImg2 = getNameArray('pa-%04d-1.jpg', 721:2:744);
        fileLoc.resize = 3/4;
        fileLoc.depthImg1 = getNameArray('result/all/demo-seq-721/plain/lagrangian/%03d_depthL.png', 1:12);
        fileLoc.depthImg2 = getNameArray('result/all/demo-seq-721/plain/lagrangian/%03d_depthR.png', 1:12);
        fileLoc.maxDepth = 10;

        testParam.maxDepth = 10;
    elseif testID == 14 || strcmp(testID, 'demo-seq-1683'),
        fileLoc.testName = 'demo-seq-1683';
        fileLoc.testID = 14;
        fileLoc.nImg = 12;
        fileLoc.inpPrefix = 'data/demo-seq-1683/';
        fileLoc.inpType = 'separate';
        fileLoc.inpImg1 = getNameArray('pa-%04d-0.jpg', 1683:2:1706);
        fileLoc.inpImg2 = getNameArray('pa-%04d-1.jpg', 1683:2:1706);
        fileLoc.resize = 3/4;
        fileLoc.depthImg1 = getNameArray('result/all/demo-seq-1683/plain/lagrangian/%03d_depthL.png', 1:12);
        fileLoc.depthImg2 = getNameArray('result/all/demo-seq-1683/plain/lagrangian/%03d_depthR.png', 1:12);
        fileLoc.maxDepth = 10;

        testParam.maxDepth = 10;
    elseif testID == 15 || strcmp(testID, 'bbb-seq-317p5'),
        fileLoc.testName = 'bbb-seq-317p5';
        fileLoc.testID = 15;
        fileLoc.nImg = 12;
        fileLoc.inpPrefix = 'data/bbb-seq-317p5/';
        fileLoc.inpType = 'separate';
        fileLoc.inpImg1 = getNameArray('bbb-317p5_%03d-0.jpg', 1:2:24);
        fileLoc.inpImg2 = getNameArray('bbb-317p5_%03d-1.jpg', 1:2:24);
        fileLoc.resize = 1;
        fileLoc.depthImg1 = getNameArray('result/all/bbb-seq-317p5/full/lagrangian/%03d_depthL.png', 1:12);
        fileLoc.depthImg2 = getNameArray('result/all/bbb-seq-317p5/full/lagrangian/%03d_depthR.png', 1:12);
        fileLoc.maxDepth = 120;

        testParam.maxDepth = 120;
    elseif testID == 16 || strcmp(testID, 'avenger-177-001'),
        fileLoc.testName = 'avenger-177-001';
        fileLoc.testID = 16;
        fileLoc.nImg = 1;
        fileLoc.inpPrefix = 'data/';
        fileLoc.inpType = 'combinedInv';
        fileLoc.inpImg = {'avenger-177_001.jpg'};
        fileLoc.resize = 1;
        fileLoc.depthImg1 = getNameArray('result/all/avenger-177-001/full-v2/lagrangian/%03d_depthL.png', 1:1);
        fileLoc.depthImg2 = getNameArray('result/all/avenger-177-001/full-v2/lagrangian/%03d_depthR.png', 1:1);
        fileLoc.maxDepth = 20;

        testParam.maxDepth = 20;
    elseif testID == 17 || strcmp(testID, 'avenger-492-001'),
        fileLoc.testName = 'avenger-492-001';
        fileLoc.testID = 17;
        fileLoc.nImg = 1;
        fileLoc.inpPrefix = 'data/';
        fileLoc.inpType = 'combinedInv';
        fileLoc.inpImg = {'avenger-492_001.jpg'};
        fileLoc.resize = 1;
        fileLoc.depthImg1 = getNameArray('result/all/avenger-492-001/full-v2/lagrangian/%03d_depthL.png', 1:1);
        fileLoc.depthImg2 = getNameArray('result/all/avenger-492-001/full-v2/lagrangian/%03d_depthR.png', 1:1);
        fileLoc.maxDepth = 20;

        testParam.maxDepth = 20;
    elseif testID == 18 || strcmp(testID, 'avenger-665-004'),
        fileLoc.testName = 'avenger-665-004';
        fileLoc.testID = 18;
        fileLoc.nImg = 1;
        fileLoc.inpPrefix = 'data/';
        fileLoc.inpType = 'combinedInv';
        fileLoc.inpImg = {'avenger-665_004.jpg'};
        fileLoc.resize = 1;
        fileLoc.depthImg1 = getNameArray('result/all/avenger-665-004/full-v2/lagrangian/%03d_depthL.png', 1:1);
        fileLoc.depthImg2 = getNameArray('result/all/avenger-665-004/full-v2/lagrangian/%03d_depthR.png', 1:1);
        fileLoc.maxDepth = 20;

        testParam.maxDepth = 20;
    elseif testID == 19 || strcmp(testID, 'despicable-29-001'),
        fileLoc.testName = 'despicable-29-001';
        fileLoc.testID = 19;
        fileLoc.nImg = 1;
        fileLoc.inpPrefix = 'data/';
        fileLoc.inpType = 'combinedInv';
        fileLoc.inpImg = {'despicable-29_001.jpg'};
        fileLoc.resize = 1;
        fileLoc.depthImg1 = getNameArray('result/all/despicable-29-001/full-v2/lagrangian/%03d_depthL.png', 1:1);
        fileLoc.depthImg2 = getNameArray('result/all/despicable-29-001/full-v2/lagrangian/%03d_depthR.png', 1:1);
        fileLoc.maxDepth = 20;

        testParam.maxDepth = 20;
    elseif testID == 20 || strcmp(testID, 'despicable-124-018'),
        fileLoc.testName = 'despicable-124-018';
        fileLoc.testID = 20;
        fileLoc.nImg = 1;
        fileLoc.inpPrefix = 'data/';
        fileLoc.inpType = 'combinedInv';
        fileLoc.inpImg = {'despicable-124_018.jpg'};
        fileLoc.resize = 1;
        fileLoc.depthImg1 = getNameArray('result/all/despicable-124-018/full-v2/lagrangian/%03d_depthL.png', 1:1);
        fileLoc.depthImg2 = getNameArray('result/all/despicable-124-018/full-v2/lagrangian/%03d_depthR.png', 1:1);
        fileLoc.maxDepth = 32;

        testParam.maxDepth = 20;
    elseif testID == 21 || strcmp(testID, 'gravity-23-001'),
        fileLoc.testName = 'gravity-23-001';
        fileLoc.testID = 21;
        fileLoc.nImg = 1;
        fileLoc.inpPrefix = 'data/';
        fileLoc.inpType = 'combinedInv';
        fileLoc.inpImg = {'gravity-23_001.jpg'};
        fileLoc.resize = 1;
        fileLoc.depthImg1 = getNameArray('result/all/gravity-23-001/full-v2/lagrangian/%03d_depthL.png', 1:1);
        fileLoc.depthImg2 = getNameArray('result/all/gravity-23-001/full-v2/lagrangian/%03d_depthR.png', 1:1);
        fileLoc.maxDepth = 20;

        testParam.maxDepth = 20;
    elseif testID == 22 || strcmp(testID, 'gravity-57-001'),
        fileLoc.testName = 'gravity-57-001';
        fileLoc.testID = 22;
        fileLoc.nImg = 1;
        fileLoc.inpPrefix = 'data/';
        fileLoc.inpType = 'combinedInv';
        fileLoc.inpImg = {'gravity-57_001.jpg'};
        fileLoc.resize = 1;
        fileLoc.depthImg1 = getNameArray('result/all/gravity-57-001/full-v2/lagrangian/%03d_depthL.png', 1:1);
        fileLoc.depthImg2 = getNameArray('result/all/gravity-57-001/full-v2/lagrangian/%03d_depthR.png', 1:1);
        fileLoc.maxDepth = 20;

        testParam.maxDepth = 20;
    elseif testID == 23 || strcmp(testID, 'pirate-112-022'),
        fileLoc.testName = 'pirate-112-022';
        fileLoc.testID = 23;
        fileLoc.nImg = 1;
        fileLoc.inpPrefix = 'data/';
        fileLoc.inpType = 'combinedInv';
        fileLoc.inpImg = {'pirate-112_022.jpg'};
        fileLoc.resize = 1;
        fileLoc.depthImg1 = getNameArray('result/all/pirate-112-022-001/full-v2/lagrangian/%03d_depthL.png', 1:1);
        fileLoc.depthImg2 = getNameArray('result/all/pirate-112-022-001/full-v2/lagrangian/%03d_depthR.png', 1:1);
        fileLoc.maxDepth = 20;

        testParam.maxDepth = 20;
    elseif testID == 24 || strcmp(testID, 'spiderman-426-024'),
        fileLoc.testName = 'spiderman-426-024';
        fileLoc.testID = 24;
        fileLoc.nImg = 1;
        fileLoc.inpPrefix = 'data/';
        fileLoc.inpType = 'combinedInv';
        fileLoc.inpImg = {'spiderman-426_024.jpg'};
        fileLoc.resize = 1;
        fileLoc.depthImg1 = getNameArray('result/all/spiderman-426-024/full-v2/lagrangian/%03d_depthL.png', 1:1);
        fileLoc.depthImg2 = getNameArray('result/all/spiderman-426-024/full-v2/lagrangian/%03d_depthR.png', 1:1);
        fileLoc.maxDepth = 20;

        testParam.maxDepth = 20;

    elseif testID == 25 || strcmp(testID, 'avenger-seq-18'),
        fileLoc.testName = 'avenger-seq-18';
        fileLoc.testID = 25;
        fileLoc.nImg = 40;
        fileLoc.inpPrefix = 'data/avenger-seq-18/';
        fileLoc.inpType = 'combinedInv';
        fileLoc.inpImg = getNameArray('avenger-18_%03d.jpg', 26:2:105);
        fileLoc.resize = 1;
        fileLoc.depthImg1 = getNameArray('result/all/avenger-seq-18/full-v2/lagrangian/%03d_depthL.png', 1:40);
        fileLoc.depthImg2 = getNameArray('result/all/avenger-seq-18/full-v2/lagrangian/%03d_depthR.png', 1:40);
        fileLoc.maxDepth = 20;

        testParam.maxDepth = 20;
    else
        error('Test case not defined');
    end
end

function nameArr = getNameArray(formatStr, inp)
    nameArr = cell(size(inp));
    for ind = 1:numel(inp);
        nameArr{ind} = sprintf(formatStr, inp(ind));
    end
end