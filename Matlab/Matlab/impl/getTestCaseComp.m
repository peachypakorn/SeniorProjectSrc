function [fileLoc, testParam] = getTestCaseComp(testID)
    fileLoc = struct();
    testParam = struct();
    if testID == 1 || strcmp(testID, 'bbb-basic'),
        fileLoc.testName = 'bbb-basic';
        fileLoc.testID = 1;
        fileLoc.nImg = 1;
        fileLoc.inpPrefix = '';
        fileLoc.inpType = 'separate';
        fileLoc.inpImg1 = {'frame1.jpg'};
        fileLoc.inpImg2 = {'frame2.jpg'};
        fileLoc.resize = 1;
        fileLoc.maxDepth = 18;

        testParam.maxDepth = 18;
    elseif testID == 2 || strcmp(testID, 'avenger-177-001'),
        fileLoc.testName = 'avenger-177-001';
        fileLoc.testID = 2;
        fileLoc.nImg = 1;
        fileLoc.inpPrefix = '';
        fileLoc.inpType = 'combinedInv';
        fileLoc.inpImg = {'avenger-177_001.jpg'};
        fileLoc.resize = 1;
        fileLoc.maxDepth = 20;

        testParam.maxDepth = 20;
    elseif testID == 3 || strcmp(testID, 'avenger-492-001'),
        fileLoc.testName = 'avenger-492-001';
        fileLoc.testID = 3;
        fileLoc.nImg = 1;
        fileLoc.inpPrefix = '';
        fileLoc.inpType = 'combinedInv';
        fileLoc.inpImg = {'avenger-492_001.jpg'};
        fileLoc.resize = 1;
        fileLoc.maxDepth = 20;

        testParam.maxDepth = 20;
    elseif testID == 4 || strcmp(testID, 'avenger-665-004'),
        fileLoc.testName = 'avenger-665-004';
        fileLoc.testID = 4;
        fileLoc.nImg = 1;
        fileLoc.inpPrefix = '';
        fileLoc.inpType = 'combinedInv';
        fileLoc.inpImg = {'avenger-665_004.jpg'};
        fileLoc.resize = 1;
        fileLoc.maxDepth = 20;

        testParam.maxDepth = 20;
    elseif testID == 5 || strcmp(testID, 'despicable-29-001'),
        fileLoc.testName = 'despicable-29-001';
        fileLoc.testID = 5;
        fileLoc.nImg = 1;
        fileLoc.inpPrefix = '';
        fileLoc.inpType = 'combinedInv';
        fileLoc.inpImg = {'despicable-29_001.jpg'};
        fileLoc.resize = 1;
        fileLoc.maxDepth = 20;

        testParam.maxDepth = 20;
    elseif testID == 6 || strcmp(testID, 'despicable-124-018'),
        fileLoc.testName = 'despicable-124-018';
        fileLoc.testID = 6;
        fileLoc.nImg = 1;
        fileLoc.inpPrefix = '';
        fileLoc.inpType = 'combinedInv';
        fileLoc.inpImg = {'despicable-124_018.jpg'};
        fileLoc.resize = 1;
        fileLoc.maxDepth = 32;

        testParam.maxDepth = 32;
    elseif testID == 7 || strcmp(testID, 'gravity-23-001'),
        fileLoc.testName = 'gravity-23-001';
        fileLoc.testID = 7;
        fileLoc.nImg = 1;
        fileLoc.inpPrefix = '';
        fileLoc.inpType = 'combinedInv';
        fileLoc.inpImg = {'gravity-23_001.jpg'};
        fileLoc.resize = 1;
        fileLoc.maxDepth = 20;

        testParam.maxDepth = 20;
    elseif testID == 8 || strcmp(testID, 'gravity-57-001'),
        fileLoc.testName = 'gravity-57-001';
        fileLoc.testID = 8;
        fileLoc.nImg = 1;
        fileLoc.inpPrefix = '';
        fileLoc.inpType = 'combinedInv';
        fileLoc.inpImg = {'gravity-57_001.jpg'};
        fileLoc.resize = 1;
        fileLoc.maxDepth = 20;

        testParam.maxDepth = 20;
    elseif testID == 9 || strcmp(testID, 'pirate-112-022'),
        fileLoc.testName = 'pirate-112-022';
        fileLoc.testID = 9;
        fileLoc.nImg = 1;
        fileLoc.inpPrefix = '';
        fileLoc.inpType = 'combinedInv';
        fileLoc.inpImg = {'pirate-112_022.jpg'};
        fileLoc.resize = 1;
        fileLoc.maxDepth = 20;

        testParam.maxDepth = 20;
    elseif testID == 10 || strcmp(testID, 'spiderman-426-024'),
        fileLoc.testName = 'spiderman-426-024';
        fileLoc.testID = 10;
        fileLoc.nImg = 1;
        fileLoc.inpPrefix = '';
        fileLoc.inpType = 'combinedInv';
        fileLoc.inpImg = {'spiderman-426_024.jpg'};
        fileLoc.resize = 1;
        fileLoc.maxDepth = 20;

        testParam.maxDepth = 20;
    % below: seq
    elseif testID == 11 || strcmp(testID, 'dracular-seq-215'),
        fileLoc.testName = 'dracular-seq-215';
        fileLoc.testID = 11;
        fileLoc.nImg = 12;
        fileLoc.inpPrefix = 'dracular-seq-215/';
        fileLoc.inpType = 'separate';
        fileLoc.inpImg1 = getNameArray('ud-%04d-0.jpg', 215:2:238);
        fileLoc.inpImg2 = getNameArray('ud-%04d-1.jpg', 215:2:238);
        fileLoc.resize = 1;
        fileLoc.maxDepth = 18;

        testParam.maxDepth = 18;
    elseif testID == 12 || strcmp(testID, 'dracular-seq-769'),
        fileLoc.testName = 'dracular-seq-769';
        fileLoc.testID = 12;
        fileLoc.nImg = 12;
        fileLoc.inpPrefix = 'dracular-seq-769/';
        fileLoc.inpType = 'separate';
        fileLoc.inpImg1 = getNameArray('ud-%04d-0.jpg', 769:2:792);
        fileLoc.inpImg2 = getNameArray('ud-%04d-1.jpg', 769:2:792);
        fileLoc.resize = 1;
        fileLoc.maxDepth = 24;

        testParam.maxDepth = 24;
    elseif testID == 13 || strcmp(testID, 'dracular-seq-969'),
        fileLoc.testName = 'dracular-seq-969';
        fileLoc.testID = 13;
        fileLoc.nImg = 12;
        fileLoc.inpPrefix = 'dracular-seq-969/';
        fileLoc.inpType = 'separate';
        fileLoc.inpImg1 = getNameArray('ud-%04d-0.jpg', 969:2:992);
        fileLoc.inpImg2 = getNameArray('ud-%04d-1.jpg', 969:2:992);
        fileLoc.resize = 1;
        fileLoc.maxDepth = 60;

        testParam.maxDepth = 60;
    elseif testID == 14 || strcmp(testID, 'airport-seq-110'),
        fileLoc.testName = 'airport-seq-110';
        fileLoc.testID = 14;
        fileLoc.nImg = 12;
        fileLoc.inpPrefix = 'airport-seq-110/airport-110-';
        fileLoc.inpType = 'combined';
        fileLoc.inpImg = getNameArray('ud-pa_%03d.jpg', 175:2:198);
        fileLoc.resize = 1;
        fileLoc.maxDepth = 27;

        testParam.maxDepth = 27;
    elseif testID == 15 || strcmp(testID, 'airport-seq-243'),
        fileLoc.testName = 'airport-seq-243';
        fileLoc.testID = 15;
        fileLoc.nImg = 12;
        fileLoc.inpPrefix = 'airport-seq-243/airport-243-';
        fileLoc.inpType = 'combined';
        fileLoc.inpImg = getNameArray('ud-pa_%03d.jpg', 21:2:44);
        fileLoc.resize = 1;
        fileLoc.maxDepth = 27;

        testParam.maxDepth = 27;
    elseif testID == 16 || strcmp(testID, 'skydiving-seq-40'),
        fileLoc.testName = 'skydiving-seq-40';
        fileLoc.testID = 16;
        fileLoc.nImg = 12;
        fileLoc.inpPrefix = 'skydiving-seq-40/skydiving-40-';
        fileLoc.inpType = 'combined';
        fileLoc.inpImg = getNameArray('ud-pa_%03d.jpg', 31:2:54);
        fileLoc.resize = 1;
        fileLoc.maxDepth = 27;

        testParam.maxDepth = 27;

    elseif testID == 17 || strcmp(testID, 'heidelberg-seq-338'),
        fileLoc.testName = 'heidelberg-seq-338';
        fileLoc.testID = 17;
        fileLoc.nImg = 12;
        fileLoc.inpPrefix = 'heidelberg-seq-338/heidelberg-338-';
        fileLoc.inpType = 'combined';
        fileLoc.inpImg = getNameArray('ud_%03d.jpg', 1:2:24);
        fileLoc.resize = 1;
        fileLoc.maxDepth = 13;

        testParam.maxDepth = 13;
    elseif testID == 18 || strcmp(testID, 'demo-seq-579'),
        fileLoc.testName = 'demo-seq-579';
        fileLoc.testID = 18;
        fileLoc.nImg = 12;
        fileLoc.inpPrefix = 'demo-seq-579/';
        fileLoc.inpType = 'separate';
        fileLoc.inpImg1 = getNameArray('%04d-0.jpg', 579:2:602);
        fileLoc.inpImg2 = getNameArray('%04d-1.jpg', 579:2:602);
        fileLoc.resize = 1;
        fileLoc.maxDepth = 13;

        testParam.maxDepth = 13;
    elseif testID == 19 || strcmp(testID, 'demo-seq-721'),
        fileLoc.testName = 'demo-seq-721';
        fileLoc.testID = 19;
        fileLoc.nImg = 12;
        fileLoc.inpPrefix = 'demo-seq-721/';
        fileLoc.inpType = 'separate';
        fileLoc.inpImg1 = getNameArray('pa-%04d-0.jpg', 721:2:744);
        fileLoc.inpImg2 = getNameArray('pa-%04d-1.jpg', 721:2:744);
        fileLoc.resize = 1;
        fileLoc.maxDepth = 13;

        testParam.maxDepth = 13;
    elseif testID == 20 || strcmp(testID, 'demo-seq-1683'),
        fileLoc.testName = 'demo-seq-1683';
        fileLoc.testID = 20;
        fileLoc.nImg = 12;
        fileLoc.inpPrefix = 'demo-seq-1683/';
        fileLoc.inpType = 'separate';
        fileLoc.inpImg1 = getNameArray('pa-%04d-0.jpg', 1683:2:1706);
        fileLoc.inpImg2 = getNameArray('pa-%04d-1.jpg', 1683:2:1706);
        fileLoc.resize = 1;
        fileLoc.maxDepth = 13;

        testParam.maxDepth = 13;
    %%% below: slow
    elseif testID == 21 || strcmp(testID, 'bbb-seq-317p5'),
        fileLoc.testName = 'bbb-seq-317p5';
        fileLoc.testID = 21;
        fileLoc.nImg = 12;
        fileLoc.inpPrefix = 'bbb-seq-317p5/';
        fileLoc.inpType = 'separate';
        fileLoc.inpImg1 = getNameArray('bbb-317p5_%03d-0.jpg', 1:2:24);
        fileLoc.inpImg2 = getNameArray('bbb-317p5_%03d-1.jpg', 1:2:24);
        fileLoc.resize = 1;
        fileLoc.maxDepth = 72;

        testParam.maxDepth = 72;

    elseif testID == 22 || strcmp(testID, 'avenger-seq-18'),
        fileLoc.testName = 'avenger-seq-18';
        fileLoc.testID = 22;
        fileLoc.nImg = 40;
        fileLoc.inpPrefix = 'avenger-seq-18/';
        fileLoc.inpType = 'combinedInv';
        fileLoc.inpImg = getNameArray('avenger-18_%03d.jpg', 26:2:105);
        fileLoc.resize = 1;
        fileLoc.maxDepth = 20;

        testParam.maxDepth = 20;

    elseif testID == 23 || strcmp(testID, 'art'),
        fileLoc.testName = 'art';
        fileLoc.testID = 23;
        fileLoc.nImg = 1;
        fileLoc.inpPrefix = '';
        fileLoc.inpType = 'separate';
        fileLoc.inpImg1 = {'art_1.png'};
        fileLoc.inpImg2 = {'art_2.png'};
        fileLoc.gtDepth = {'art_d.png'};
        fileLoc.resize = 1;
        fileLoc.maxDepth = 30;
        fileLoc.gtMul = -1/4;
        fileLoc.gtAdd = 0;
        fileLoc.cut = -30;

        testParam.maxDepth = 30;

    elseif testID == 24 || strcmp(testID, 'doll'),
        fileLoc.testName = 'doll';
        fileLoc.testID = 24;
        fileLoc.nImg = 1;
        fileLoc.inpPrefix = '';
        fileLoc.inpType = 'separate';
        fileLoc.inpImg1 = {'doll_1.png'};
        fileLoc.inpImg2 = {'doll_2.png'};
        fileLoc.gtDepth = {'doll_d.png'};
        fileLoc.resize = 1;
        fileLoc.maxDepth = 30;
        fileLoc.gtMul = -1/4;
        fileLoc.gtAdd = 0;
        fileLoc.cut = -30;

        testParam.maxDepth = 30;

    elseif testID == 25 || strcmp(testID, 'aloe'),
        fileLoc.testName = 'aloe';
        fileLoc.testID = 25;
        fileLoc.nImg = 1;
        fileLoc.inpPrefix = '';
        fileLoc.inpType = 'separate';
        fileLoc.inpImg1 = {'aloe_1.png'};
        fileLoc.inpImg2 = {'aloe_2.png'};
        fileLoc.resize = 1;
        fileLoc.maxDepth = 45;
        fileLoc.gtMul = -1/4;
        fileLoc.gtAdd = 0;
        fileLoc.cut = -30;

        testParam.maxDepth = 45;

    elseif testID == 26 || strcmp(testID, 'ed-301_041'),
        fileLoc.testName = 'ed-301_041';
        fileLoc.testID = 26;
        fileLoc.nImg = 1;
        fileLoc.inpPrefix = '';
        fileLoc.inpType = 'combined';
        fileLoc.inpImg = {'ed-301_041.jpg'};
        fileLoc.resize = 1;
        fileLoc.maxDepth = 32;

        testParam.maxDepth = 32;

    elseif testID == 27 || strcmp(testID, 'ed-321_025'),
        fileLoc.testName = 'ed-321_025';
        fileLoc.testID = 27;
        fileLoc.nImg = 1;
        fileLoc.inpPrefix = '';
        fileLoc.inpType = 'combined';
        fileLoc.inpImg = {'ed-321_025.jpg'};
        fileLoc.resize = 1;
        fileLoc.maxDepth = 32;

        testParam.maxDepth = 32;
    elseif testID == 28 || strcmp(testID, 'ed-443_026'),
        fileLoc.testName = 'ed-443_026';
        fileLoc.testID = 28;
        fileLoc.nImg = 1;
        fileLoc.inpPrefix = '';
        fileLoc.inpType = 'combined';
        fileLoc.inpImg = {'ed-443_026.jpg'};
        fileLoc.resize = 1;
        fileLoc.maxDepth = 32;

        testParam.maxDepth = 32;
    elseif testID == 29 || strcmp(testID, 'ed-seq-118'),
        fileLoc.testName = 'ed-seq-118';
        fileLoc.testID = 29;
        fileLoc.nImg = 12;
        fileLoc.inpPrefix = 'ed-seq-118/';
        fileLoc.inpType = 'combined';
        fileLoc.inpImg = getNameArray('ed-118_%03d.jpg', 1:2:24);
        fileLoc.resize = 1;
        fileLoc.maxDepth = 32;

        testParam.maxDepth = 32;

    elseif testID == 30 || strcmp(testID, 'ed-48_001'),
        fileLoc.testName = 'ed-48_001';
        fileLoc.testID = 30;
        fileLoc.nImg = 1;
        fileLoc.inpPrefix = '';
        fileLoc.inpType = 'combined';
        fileLoc.inpImg = {'ed-48_001.jpg'};
        fileLoc.resize = 1;
        fileLoc.maxDepth = 32;

        testParam.maxDepth = 32;
    elseif testID == 31 || strcmp(testID, 'ed-506_014'),
        fileLoc.testName = 'ed-506_014';
        fileLoc.testID = 31;
        fileLoc.nImg = 1;
        fileLoc.inpPrefix = '';
        fileLoc.inpType = 'combined';
        fileLoc.inpImg = {'ed-506_014.jpg'};
        fileLoc.resize = 1;
        fileLoc.maxDepth = 48;

        testParam.maxDepth = 48;
    elseif testID == 32 || strcmp(testID, 'skullrock_1232'),
        fileLoc.testName = 'skullrock_1232';
        fileLoc.testID = 32;
        fileLoc.nImg = 1;
        fileLoc.inpPrefix = '';
        fileLoc.inpType = 'topDown';
        fileLoc.inpImg = {'skullrock_1232.jpg'};
        fileLoc.resize = 1;
        fileLoc.maxDepth = 76;
        fileLoc.cut = -12;

        testParam.maxDepth = 76;
    elseif testID == 33 || strcmp(testID, 'skullrock_1344'),
        fileLoc.testName = 'skullrock_1344';
        fileLoc.testID = 33;
        fileLoc.nImg = 1;
        fileLoc.inpPrefix = '';
        fileLoc.inpType = 'topDown';
        fileLoc.inpImg = {'skullrock_1344.jpg'};
        fileLoc.resize = 1;
        fileLoc.maxDepth = 168;

        testParam.maxDepth = 168;
    elseif testID == 34 || strcmp(testID, 'skullrock_1478'),
        fileLoc.testName = 'skullrock_1478';
        fileLoc.testID = 34;
        fileLoc.nImg = 1;
        fileLoc.inpPrefix = '';
        fileLoc.inpType = 'topDown';
        fileLoc.inpImg = {'skullrock_1478.jpg'};
        fileLoc.resize = 1;
        fileLoc.maxDepth = 48;
        fileLoc.cut = 28;

        testParam.maxDepth = 48;
    elseif testID == 35 || strcmp(testID, 'dracular_1125'),
        fileLoc.testName = 'dracular_1125';
        fileLoc.testID = 35;
        fileLoc.nImg = 1;
        fileLoc.inpPrefix = '';
        fileLoc.inpType = 'separate';
        fileLoc.inpImg1 = {'dracular_1125-0.jpg'};
        fileLoc.inpImg2 = {'dracular_1125-1.jpg'};
        fileLoc.resize = 1;
        fileLoc.maxDepth = 44;

        testParam.maxDepth = 44;
    elseif testID == 36 || strcmp(testID, 'dracular_1225'),
        fileLoc.testName = 'dracular_1225';
        fileLoc.testID = 36;
        fileLoc.nImg = 1;
        fileLoc.inpPrefix = '';
        fileLoc.inpType = 'separate';
        fileLoc.inpImg1 = {'dracular_1225-0.jpg'};
        fileLoc.inpImg2 = {'dracular_1225-1.jpg'};
        fileLoc.resize = 1;
        fileLoc.maxDepth = 44;

        testParam.maxDepth = 44;

    elseif testID == 37 || strcmp(testID, 'demo_1003'),
        fileLoc.testName = 'demo_1003';
        fileLoc.testID = 37;
        fileLoc.nImg = 1;
        fileLoc.inpPrefix = '';
        fileLoc.inpType = 'separate';
        fileLoc.inpImg1 = {'demo_1003-0.jpg'};
        fileLoc.inpImg2 = {'demo_1003-1.jpg'};
        fileLoc.resize = 1;
        fileLoc.maxDepth = 20;

        testParam.maxDepth = 20;

    elseif testID == 38 || strcmp(testID, 'demo_1081'),
        fileLoc.testName = 'demo_1081';
        fileLoc.testID = 38;
        fileLoc.nImg = 1;
        fileLoc.inpPrefix = '';
        fileLoc.inpType = 'separate';
        fileLoc.inpImg1 = {'demo_1081-0.jpg'};
        fileLoc.inpImg2 = {'demo_1081-1.jpg'};
        fileLoc.resize = 1;
        fileLoc.maxDepth = 20;

        testParam.maxDepth = 20;

    elseif testID == 39 || strcmp(testID, 'bbb-81_050'),
        fileLoc.testName = 'bbb-81_050';
        fileLoc.testID = 39;
        fileLoc.nImg = 1;
        fileLoc.inpPrefix = '';
        fileLoc.inpType = 'separate';
        fileLoc.inpImg1 = {'bbb-R-81_050.jpg'};
        fileLoc.inpImg2 = {'bbb-L-81_050.jpg'};
        fileLoc.resize = 1;
        fileLoc.maxDepth = 72;

        testParam.maxDepth = 72;

    elseif testID == 40 || strcmp(testID, 'bbb-265_006'),
        fileLoc.testName = 'bbb-265_006';
        fileLoc.testID = 40;
        fileLoc.nImg = 1;
        fileLoc.inpPrefix = '';
        fileLoc.inpType = 'separate';
        fileLoc.inpImg1 = {'bbb-R-265_006.jpg'};
        fileLoc.inpImg2 = {'bbb-L-265_006.jpg'};
        fileLoc.resize = 1;
        fileLoc.maxDepth = 32;

        testParam.maxDepth = 32;

    elseif testID == 41 || strcmp(testID, 'bbb-252p7'),
        fileLoc.testName = 'bbb-252p7';
        fileLoc.testID = 41;
        fileLoc.nImg = 1;
        fileLoc.inpPrefix = '';
        fileLoc.inpType = 'separate';
        fileLoc.inpImg1 = {'bbb_252p7_right.jpg'};
        fileLoc.inpImg2 = {'bbb_252p7_left.jpg'};
        fileLoc.resize = 1;
        fileLoc.maxDepth = 60;

        testParam.maxDepth = 60;

    elseif testID == 42 || strcmp(testID, 'bbb-106'),
        fileLoc.testName = 'bbb-106';
        fileLoc.testID = 42;
        fileLoc.nImg = 1;
        fileLoc.inpPrefix = '';
        fileLoc.inpType = 'separate';
        fileLoc.inpImg1 = {'bbb_106_right.jpg'};
        fileLoc.inpImg2 = {'bbb_106_left.jpg'};
        fileLoc.resize = 1;
        fileLoc.maxDepth = 60;

        testParam.maxDepth = 60;

    elseif testID == 43 || strcmp(testID, 'transparency'),
        fileLoc.testName = 'transparency';
        fileLoc.testID = 43;
        fileLoc.nImg = 1;
        fileLoc.inpPrefix = '';
        fileLoc.inpType = 'separate';
        fileLoc.inpImg1 = {'transparency_right.jpg'};
        fileLoc.inpImg2 = {'transparency_left.jpg'};
        fileLoc.resize = 1;
        fileLoc.maxDepth = 40;

        testParam.maxDepth = 40;

    elseif testID == 44 || strcmp(testID, 'teapot'),
        fileLoc.testName = 'teapot';
        fileLoc.testID = 44;
        fileLoc.nImg = 1;
        fileLoc.inpPrefix = '';
        fileLoc.inpType = 'separate';
        fileLoc.inpImg1 = {'teapot_16.jpg'};
        fileLoc.inpImg2 = {'teapot_27.jpg'};
        fileLoc.resize = 1;
        fileLoc.maxDepth = 24;

        testParam.maxDepth = 24;

    elseif testID == 45 || strcmp(testID, 'ed-seq2s-118'),
        fileLoc.testName = 'ed-seq2s-118';
        fileLoc.testID = 45;
        fileLoc.nImg = 49;
        fileLoc.inpPrefix = 'ed-seq-118/';
        fileLoc.inpType = 'combined';
        fileLoc.inpImg = getNameArray('ed-118_%03d.jpg', 1:1:49);
        fileLoc.resize = 0.5;
        fileLoc.maxDepth = 16;

        testParam.maxDepth = 16;

    elseif testID == 46 || strcmp(testID, 'bbb-seq2s-265'),
        fileLoc.testName = 'bbb-seq2s-265';
        fileLoc.testID = 46;
        fileLoc.nImg = 49;
        fileLoc.inpPrefix = 'bbb-seq-265/';
        fileLoc.inpType = 'separate';
        fileLoc.inpImg1 = getNameArray('bbb-265-R-265_%03d.jpg', 1:1:49);
        fileLoc.inpImg2 = getNameArray('bbb-265-L-265_%03d.jpg', 1:1:49);
        fileLoc.resize = 0.5;
        fileLoc.maxDepth = 16;

        testParam.maxDepth = 16;

    elseif testID == 47 || strcmp(testID, 'bbb-265-ld'),
        fileLoc.testName = 'bbb-265-ld';
        fileLoc.testID = 47;
        fileLoc.nImg = 1;
        fileLoc.inpPrefix = '';
        fileLoc.inpType = 'separate';
        fileLoc.inpImg1 = {'bbb-265-R-265_001.jpg'};
        fileLoc.inpImg2 = {'bbb-265-L-265_001.jpg'};
        fileLoc.resize = 1;
        fileLoc.maxDepth = 32;

        testParam.maxDepth = 32;

    elseif testID == 48 || strcmp(testID, 'ed-118-ld'),
        fileLoc.testName = 'ed-118-ld';
        fileLoc.testID = 48;
        fileLoc.nImg = 1;
        fileLoc.inpPrefix = '';
        fileLoc.inpType = 'combined';
        fileLoc.inpImg = getNameArray('ed-118_%03d.jpg', 1:1:1);
        fileLoc.resize = 1;
        fileLoc.maxDepth = 32;

        testParam.maxDepth = 32;

    elseif testID == 49 || strcmp(testID, 'demo-26_001'),
        fileLoc.testName = 'demo-26_001';
        fileLoc.testID = 49;
        fileLoc.nImg = 1;
        fileLoc.inpPrefix = '';
        fileLoc.inpType = 'topDown';
        fileLoc.inpImg = {'demo-26_001.jpg'};
        fileLoc.resize = 1;
        fileLoc.maxDepth = 32;

        testParam.maxDepth = 32;
    elseif testID == 50 || strcmp(testID, 'demo-seq-26'),
        fileLoc.testName = 'demo-seq-26';
        fileLoc.testID = 50;
        fileLoc.nImg = 48;
        fileLoc.inpPrefix = 'demo2-seq-26/';
        fileLoc.inpType = 'topDown';
        fileLoc.inpImg = getNameArray('demo-26_%03d.jpg', 1:1:48);;
        fileLoc.resize = 0.5;
        fileLoc.maxDepth = 16;

        testParam.maxDepth = 16;
    elseif testID == 51 || strcmp(testID, 'skullrock-seq-1478'),
        fileLoc.testName = 'skullrock-seq-1478';
        fileLoc.testID = 51;
        fileLoc.nImg = 39;
        fileLoc.inpPrefix = 'skullrock-seq-1478/';
        fileLoc.inpType = 'topDown';
        fileLoc.inpImg = getNameArray('skullrock_%03d.jpg', 1458:1:1496);
        fileLoc.resize = 0.5;
        fileLoc.maxDepth = 28;
        fileLoc.cut = 14;

        testParam.maxDepth = 28;
    elseif testID == 52 || strcmp(testID, 'skullrock-1496'),
        fileLoc.testName = 'skullrock-seq-1496';
        fileLoc.testID = 52;
        fileLoc.nImg = 1;
        fileLoc.inpPrefix = '';
        fileLoc.inpType = 'separate';
        fileLoc.inpImg1 = {'skullrock_L_1496.jpg'};
        fileLoc.inpImg2 = {'skullrock_R_1496.jpg'};
        fileLoc.resize = 1;
        fileLoc.maxDepth = 48;
        fileLoc.cut = 28;

        testParam.maxDepth = 48;

    elseif testID == 53 || strcmp(testID, 'bbb-190_131'),
        fileLoc.testName = 'bbb-190_131';
        fileLoc.testID = 53;
        fileLoc.nImg = 1;
        fileLoc.inpPrefix = '';
        fileLoc.inpType = 'separate';
        fileLoc.inpImg1 = {'bbb_R-190_131.jpg'};
        fileLoc.inpImg2 = {'bbb_L-190_131.jpg'};
        fileLoc.resize = 1;
        fileLoc.maxDepth = 32;

    elseif testID == 54 || strcmp(testID, 'ed-200_200'),
        fileLoc.testName = 'ed-200_200';
        fileLoc.testID = 54;
        fileLoc.nImg = 1;
        fileLoc.inpPrefix = '';
        fileLoc.inpType = 'combined';
        fileLoc.inpImg = {'ed-200_200.jpg'};
        fileLoc.resize = 1;
        fileLoc.maxDepth = 32;

        testParam.maxDepth = 32;

    elseif testID == 55 || strcmp(testID, 'bbb-321_084'),
        fileLoc.testName = 'bbb-321_084';
        fileLoc.testID = 55;
        fileLoc.nImg = 1;
        fileLoc.inpPrefix = '';
        fileLoc.inpType = 'separate';
        fileLoc.inpImg1 = {'bbb_R-321_084.jpg'};
        fileLoc.inpImg2 = {'bbb_L-321_084.jpg'};
        fileLoc.resize = 1;
        fileLoc.maxDepth = 48;
        
    elseif testID == 56 || strcmp(testID, 'sintel-262_001'),
        fileLoc.testName = 'sintel-262_001';
        fileLoc.testID = 56;
        fileLoc.nImg = 1;
        fileLoc.inpPrefix = '';
        fileLoc.inpType = 'combined';
        fileLoc.inpImg = {'sintel-262_001.jpg'};
        fileLoc.resize = 1;
        fileLoc.maxDepth = 24;

        testParam.maxDepth = 24;
        
    elseif testID == 57 || strcmp(testID, 'sintel-262_035'),
        fileLoc.testName = 'sintel-262_035';
        fileLoc.testID = 57;
        fileLoc.nImg = 1;
        fileLoc.inpPrefix = '';
        fileLoc.inpType = 'combined';
        fileLoc.inpImg = {'sintel-262_035.jpg'};
        fileLoc.resize = 1;
        fileLoc.maxDepth = 24;

        testParam.maxDepth = 24;
        
    elseif testID == 58 || strcmp(testID, 'sintel-262_109'),
        fileLoc.testName = 'sintel-262_109';
        fileLoc.testID = 58;
        fileLoc.nImg = 1;
        fileLoc.inpPrefix = '';
        fileLoc.inpType = 'combined';
        fileLoc.inpImg = {'sintel-262_109.jpg'};
        fileLoc.resize = 1;
        fileLoc.maxDepth = 24;

        testParam.maxDepth = 24;
        
    elseif testID == 59 || strcmp(testID, 'sintel-262_147'),
        fileLoc.testName = 'sintel-262_147';
        fileLoc.testID = 59;
        fileLoc.nImg = 1;
        fileLoc.inpPrefix = '';
        fileLoc.inpType = 'combined';
        fileLoc.inpImg = {'sintel-262_147.jpg'};
        fileLoc.resize = 1;
        fileLoc.maxDepth = 24;

        testParam.maxDepth = 24;

    elseif testID == 60 || strcmp(testID, 'reindeer'),
        fileLoc.testName = 'reindeer';
        fileLoc.testID = 60;
        fileLoc.nImg = 1;
        fileLoc.inpPrefix = '';
        fileLoc.inpType = 'separate';
        fileLoc.inpImg1 = {'reindeer_1.png'};
        fileLoc.inpImg2 = {'reindeer_2.png'};
        fileLoc.gtDepth = {'reindeer_d.png'};
        fileLoc.resize = 1;
        fileLoc.maxDepth = 30;
        fileLoc.gtMul = -1/4;
        fileLoc.gtAdd = 0;
        fileLoc.cut = -30;

        testParam.maxDepth = 30;

    elseif testID == 61 || strcmp(testID, 'testBox1'),
        fileLoc.testName = 'testBox1';
        fileLoc.testID = 61;
        fileLoc.nImg = 1;
        fileLoc.inpPrefix = '';
        fileLoc.inpType = 'separate';
        fileLoc.inpImg1 = {'testBox1-R.png'};
        fileLoc.inpImg2 = {'testBox1-M.png'};
        fileLoc.resize = 1;
        fileLoc.maxDepth = 12;

        testParam.maxDepth = 12;

    elseif testID == 62 || strcmp(testID, 'testBox6'),
        fileLoc.testName = 'testBox6';
        fileLoc.testID = 62;
        fileLoc.nImg = 1;
        fileLoc.inpPrefix = '';
        fileLoc.inpType = 'separate';
        fileLoc.inpImg1 = {'testBox6-hd-R.png'};
        fileLoc.inpImg2 = {'testBox6-hd-M.png'};
        fileLoc.gtDepth = {'testBox6-hd-Dp.png'};
        fileLoc.resize = 1;
        fileLoc.maxDepth = 16;
        fileLoc.gtMul = 0.125;
        fileLoc.gtAdd = -128;
        fileLoc.cut = 0;

        testParam.maxDepth = 16;

    elseif testID == 63 || strcmp(testID, 'synCamBlur1')
        fileLoc.testName = 'synCamBlur1';
        fileLoc.testID = 63;
        fileLoc.nImg = 1;
        fileLoc.inpPrefix = '';
        fileLoc.inpType = 'separate';
        fileLoc.inpImg1 = {'syn_blur_Cam_A1.jpg'};
        fileLoc.inpImg2 = {'syn_blur_Cam_A5.jpg'};
        fileLoc.resize = 1;
        fileLoc.maxDepth = 80;

    elseif testID == 64 || strcmp(testID, 'avenger-seq-01267'),
        fileLoc.testName = 'avenger-seq-01267';
        fileLoc.testID = 64;
        fileLoc.nImg = 30;
        fileLoc.inpPrefix = 'avenger-seq-01267/';
        fileLoc.inpType = 'combinedInv';
        fileLoc.inpImg = getNameArray('Avenger_%05d.jpg', 1267:1:1296);
        fileLoc.resize = [540 960];
        fileLoc.maxDepth = 24;

        testParam.maxDepth = 24;

    elseif testID == 65 || strcmp(testID, 'avenger-seq-07830'),
        fileLoc.testName = 'avenger-seq-07830';
        fileLoc.testID = 65;
        fileLoc.nImg = 49;
        fileLoc.inpPrefix = 'avenger-seq-07830/';
        fileLoc.inpType = 'combinedInv';
        fileLoc.inpImg = getNameArray('Avenger_%05d.jpg', 7830:1:7878);
        fileLoc.resize = [540 960];
        fileLoc.maxDepth = 24;

        testParam.maxDepth = 24;

    elseif testID == 66 || strcmp(testID, 'pirate-seq-00748'),
        fileLoc.testName = 'pirate-seq-00748';
        fileLoc.testID = 66;
        fileLoc.nImg = 27;
        fileLoc.inpPrefix = 'pirate-seq-00748/';
        fileLoc.inpType = 'combinedInv';
        fileLoc.inpImg = getNameArray('Pirate_%05d.jpg', 748:1:774);
        fileLoc.resize = [540 960];
        fileLoc.maxDepth = 24;

        testParam.maxDepth = 24;

    elseif testID == 67 || strcmp(testID, 'pirate-seq-01425'),
        fileLoc.testName = 'pirate-seq-01425';
        fileLoc.testID = 67;
        fileLoc.nImg = 37;
        fileLoc.inpPrefix = 'pirate-seq-01425/';
        fileLoc.inpType = 'combinedInv';
        fileLoc.inpImg = getNameArray('Pirate_%05d.jpg', 1425:1:1461);
        fileLoc.resize = [540 960];
        fileLoc.maxDepth = 24;

        testParam.maxDepth = 24;

    elseif testID == 68 || strcmp(testID, 'michaeljackson-seq-01943'),
        fileLoc.testName = 'michaeljackson-seq-01943';
        fileLoc.testID = 68;
        fileLoc.nImg = 49;
        fileLoc.inpPrefix = 'michaeljackson-seq-01943/';
        fileLoc.inpType = 'combinedInv';
        fileLoc.inpImg = getNameArray('MichaelJackson_%05d.jpg', 1943:1:1991);
        fileLoc.resize = [540 960];
        fileLoc.maxDepth = 24;

        testParam.maxDepth = 24;

    elseif testID == 69 || strcmp(testID, 'michaeljackson-seq-13124'),
        fileLoc.testName = 'michaeljackson-seq-13124';
        fileLoc.testID = 69;
        fileLoc.nImg = 49;
        fileLoc.inpPrefix = 'michaeljackson-seq-13124/';
        fileLoc.inpType = 'combinedInv';
        fileLoc.inpImg = getNameArray('MichaelJackson_%05d.jpg', 13124:1:13172);
        fileLoc.resize = [540 960];
        fileLoc.maxDepth = 30;

        testParam.maxDepth = 30;

    elseif testID == 70 || strcmp(testID, 'spiderman-seq-05237'),
        fileLoc.testName = 'spiderman-seq-05237';
        fileLoc.testID = 70;
        fileLoc.nImg = 41;
        fileLoc.inpPrefix = 'spiderman-seq-05237/';
        fileLoc.inpType = 'combinedInv';
        fileLoc.inpImg = getNameArray('SpiderMan_%05d.jpg', 5237:1:5277);
        fileLoc.resize = [540 960];
        fileLoc.maxDepth = 24;

        testParam.maxDepth = 24;

    elseif testID == 71 || strcmp(testID, 'spiderman-seq-12366'),
        fileLoc.testName = 'spiderman-seq-12366';
        fileLoc.testID = 71;
        fileLoc.nImg = 49;
        fileLoc.inpPrefix = 'spiderman-seq-12366/';
        fileLoc.inpType = 'combinedInv';
        fileLoc.inpImg = getNameArray('SpiderMan_%05d.jpg', 12366:1:12414);
        fileLoc.resize = [540 960];
        fileLoc.maxDepth = 24;

        testParam.maxDepth = 24;

    elseif testID == 72 || strcmp(testID, 'dracular-seq-463'),
        fileLoc.testName = 'dracular-seq-463';
        fileLoc.testID = 72;
        fileLoc.nImg = 49;
        fileLoc.inpPrefix = 'dracular-seq-463/';
        fileLoc.inpType = 'topDown';
        fileLoc.inpImg = getNameArray('dracula-0_%03d.jpg', 463:1:511);
        fileLoc.resize = 0.5;
        fileLoc.maxDepth = 24;
        fileLoc.cut = 24;

        testParam.maxDepth = 24;

    elseif testID == 73 || strcmp(testID, 'dracular-seq-1195'),
        fileLoc.testName = 'dracular-seq-1195';
        fileLoc.testID = 73;
        fileLoc.nImg = 45;
        fileLoc.inpPrefix = 'dracular-seq-1195/';
        fileLoc.inpType = 'topDown';
        fileLoc.inpImg = getNameArray('dracula-0_%03d.jpg', 1195:1:1239);
        fileLoc.resize = 0.5;
        fileLoc.maxDepth = 80;
        fileLoc.cut = -35;

        testParam.maxDepth = 80;

    elseif testID == 74 || strcmp(testID, 'demo-seq-35'),
        fileLoc.testName = 'demo-seq-35';
        fileLoc.testID = 74;
        fileLoc.nImg = 49;
        fileLoc.inpPrefix = 'demo-seq-35/';
        fileLoc.inpType = 'topDown';
        fileLoc.inpImg = getNameArray('demo-35_%03d.jpg', 10:1:58);
        fileLoc.resize = 0.5;
        fileLoc.maxDepth = 12;

        testParam.maxDepth = 12;

    elseif testID == 75 || strcmp(testID, 'bbb-seq-190'),
        fileLoc.testName = 'bbb-seq-190';
        fileLoc.testID = 75;
        fileLoc.nImg = 49;
        fileLoc.inpPrefix = 'bbb-seq-190/';
        fileLoc.inpType = 'separate';
        fileLoc.inpImg1 = getNameArray('bbb-R-190_%03d.jpg', 101:1:149);
        fileLoc.inpImg2 = getNameArray('bbb-L-190_%03d.jpg', 101:1:149);
        fileLoc.resize = 0.5;
        fileLoc.maxDepth = 24;

        testParam.maxDepth = 24;

    elseif testID == 76 || strcmp(testID, 'ed-seq-398'),
        fileLoc.testName = 'ed-seq-398';
        fileLoc.testID = 76;
        fileLoc.nImg = 49;
        fileLoc.inpPrefix = 'ed-seq-398/';
        fileLoc.inpType = 'combined';
        fileLoc.inpImg = getNameArray('ed-398_%03d.jpg', 75:1:123);
        fileLoc.resize = [540 960];
        fileLoc.maxDepth = 16;

        testParam.maxDepth = 16;

    elseif testID == 77 || strcmp(testID, 'skullrock-seq-t33'),
        fileLoc.testName = 'skullrock-seq-t33';
        fileLoc.testID = 77;
        fileLoc.nImg = 35;
        fileLoc.inpPrefix = 'skullrock-seq-t33/';
        fileLoc.inpType = 'topDown';
        fileLoc.inpImg = getNameArray('skullrock-33_%03d.jpg', 46:1:80);
        fileLoc.resize = 0.5;
        fileLoc.maxDepth = 80;
        fileLoc.cut = -20;

        testParam.maxDepth = 80;

    elseif testID == 78 || strcmp(testID, 'pirate-seq-00758'),
        fileLoc.testName = 'pirate-seq-00758';
        fileLoc.testID = 78;
        fileLoc.nImg = 39;
        fileLoc.inpPrefix = 'pirate-seq-00758/';
        fileLoc.inpType = 'combinedInv';
        fileLoc.inpImg = getNameArray('Pirate_%05d.jpg', 758:1:796);
        fileLoc.resize = [540 960];
        fileLoc.maxDepth = 24;

        testParam.maxDepth = 24;
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