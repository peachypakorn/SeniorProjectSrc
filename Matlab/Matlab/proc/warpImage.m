function warp = warpImage(imO, vxO, depthO, nPP)
    [heightO,widthO,nchannels]=size(imO);

    warp = 0;
    for iPP = 0:nPP,
        height = ceil(heightO / (2^(nPP-iPP)));
        width = ceil(widthO / (2^(nPP-iPP)));
        im = imresize(imO, [height, width]);
        vx = imresize(vxO, [height, width]) * width / widthO;
        depth = imresize(depthO, [height, width]) * width / widthO;

        [xx,yy]=meshgrid(1:width,1:height);
        data = zeros(widthO*heightO*4, 3+nchannels);
        nD = 0;
        for iR = 1:height,
            cC = 0;
            cD = depth(iR, 1);
            cV = reshape(im(iR, 1, :), [1, nchannels]);
            for iC = 1:width+1,
                if iC ~= width+1,
                    uC = iC - vx(iR, iC);
                    uD = depth(iR, iC);
                    uV = reshape(im(iR, iC, :), [1, nchannels]);
                else,
                    uC = width + 1;
                    uD = depth(iR, width);
                    uV = reshape(im(iR, width, :), [1, nchannels]);
                end

                if abs(uD - cD) <= 2 && uC ~= cC,
                    if uC > cC,
                        rC = (uC:-1:cC+1e-10)';
                    else,
                        rC = (uC:1:cC-1e-10)';
                    end
                    rD = ((rC - uC) * cD + (cC - rC) * uD) / (cC - uC);
                    rV = ((rC - uC) * cV + (cC - rC) * uV) / (cC - uC);
                    rL = length(rC);
                    data(nD+1:nD+rL, :) = [rC, iR*ones(rL, 1), -rD, rV];
                    nD = nD + rL;
                else,
                    data(nD+1, :) = [cC, iR, -cD, cV];
                    nD = nD + 1;
                end

                cC = uC;
                cD = uD;
                cV = uV;
            end
            data(nD+1, :) = [cC, iR, -cD, cV];
            nD = nD + 1;
        end
        if iPP ~= 0,
            warp = imresize(warp, [height, width]);
            data(nD+1:nD+height*width, :) = [xx(:), yy(:), -1e10*ones(height*width, 1), ...
                                             reshape(warp, [height*width, nchannels])];
            nD = nD + height * width;
        end

        data = data(1:nD, :);
        data = [round(data(:, 1:2)), data(:, 3), data];

        data = sortrows(data);
        [~, ia, ~] = unique(data(:,1:2),'rows', 'legacy');
        data = data(ia,:);

        warp = zeros([height, width, nchannels]);
        for ind=1:nchannels
            a = data(:,4);
            b = data(:,5);
            c = data(:,6+ind);
            F = TriScatteredInterp(a,b,c);

            warp(:,:,ind) = F(xx,yy);
        end
        %if sum(abs(vxO+depthO*3)) < 1,
        %    figure;
        %    imshow(warp)
        %end
    end
    %error('123')
end