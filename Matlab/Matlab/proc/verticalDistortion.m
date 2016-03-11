function [imgC1, imgC2] = verticalDistortion(img1, img2, htor, vtor)
    imgBW1 = rgb2gray(img1);
    imgBW2 = rgb2gray(img2);
    points1 = detectSURFFeatures(imgBW1);
    points2 = detectSURFFeatures(imgBW2);
    [features1, points1] = extractFeatures(imgBW1, points1);
    [features2, points2] = extractFeatures(imgBW2, points2);
    pairs = matchFeatures(features1, features2);
    pairs = pairs(1:min(200, length(pairs)), :);
    points1 = points1.Location;
    points2 = points2.Location;
    pairs = pairs(abs(points1(pairs(:, 1), 1) - points2(pairs(:, 2), 1)) < htor & ...
                  abs(points1(pairs(:, 1), 2) - points2(pairs(:, 2), 2)) < vtor, :);

    matchedPoints1 = points1(pairs(:, 1), :);
    matchedPoints2 = points2(pairs(:, 2), :);

    A = [matchedPoints2, ones(length(matchedPoints2), 1)];
    %A = [matchedPoints2, ones(length(matchedPoints2), 1), matchedPoints2(:, 1).^2 / 1000, ...
    %     matchedPoints2(:, 2).^2 / 1000,  matchedPoints2(:, 1) .* matchedPoints2(:, 2) / 1000];
    b = matchedPoints1(:, 2);
    %v = A \ b;
    %v(4:6) = v(4:6) / 1000;
    %forward = @(U, e) [U(:, 1), U(:, 1) * v(1) + U(:, 2) * v(2) + v(3) + U(:, 1).^2 * v(4) + U(:, 2).^2 * v(5) + U(:, 1).*U(:, 2)*v(6)];
    %reverse = @(X, e) [X(:, 1), (-v(2)-X(:, 1)*v(6)+sqrt((v(2)+X(:, 1)*v(6)).^2-4*v(4)*(X(:, 1) * v(1)+v(3)+ X(:, 1).^2 * v(4) - X(:, 2))))/(2*v(4))];
    %T = maketform('custom', 2, 2, forward, reverse, v);
    Ta = zeros(3);
    Ta(1, 1) = 1;
    Ta(3, 3) = 1;
    Ta(:, 2) = A \ b;
    T = maketform('affine', Ta);

    if false
        figure
        tmp = zeros(size(img1));
        tmp(:, :, 1) = imgBW1;
        tmp(:, :, 2) = imtransform(imgBW2, T, 'XData', [1 size(img2,2)], 'YData', [1 size(img1,1)]);
        tmp(:, :, 3) = tmp(:, :, 2);
        imshow(tmp)

        figure
        tmp(:, :, 2) = imgBW2;
        tmp(:, :, 3) = imgBW2;
        imshow(tmp)
    end

    imgC1 = img1;
    imgC2 = imtransform(img2, T, 'XData', [1 size(img2,2)], 'YData', [1 size(img1,1)]);
end