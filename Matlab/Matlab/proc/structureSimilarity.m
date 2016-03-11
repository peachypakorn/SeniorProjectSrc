function [res, res_sig] = structureSimilarity(img1, img2, nH) 
	if nargin < 3,
		nH = 4;
	end
	if size(img1, 3) > 1,
		img1 = img1(:, :, 1) * 0.3 + img1(:, :, 2) * 0.6 + img1(:, :, 3) * 0.1;
		img2 = img2(:, :, 1) * 0.3 + img2(:, :, 2) * 0.6 + img2(:, :, 3) * 0.1;
	end
	nR = size(img1, 1);
	nC = size(img1, 2);
	mUx = boxFilter(img1, nH);
	mUy = boxFilter(img2, nH);
	mUxy = mUx .* mUy;
	mSxx = boxFilter(img1 .^ 2, nH) - mUx .^2;
	mSyy = boxFilter(img2 .^ 2, nH) - mUy .^2;
	mSxy = boxFilter(img1 .* img2, nH) - mUxy;

	c1 = 0.01^2;
	c2 = 0.03^2;
	res = (2 * mUxy + c1) .* (2 * mSxy + c2) ./ ((mUx.^2 + mUy.^2 + c1) .* (mSxx + mSyy + c2));
	res_sig = (2 * mSxy + c2) ./ (mSxx + mSyy + c2);
end

function res = boxFilter(img, nH)
	nR = size(img, 1);
	nC = size(img, 2);
	mCnt = zeros(nR, nC);
	sx = zeros(nR, nC);
	res = zeros(nR, nC);
	for dx = -nH:nH,
		iL = max(1+dx, 1) - dx;
		iR = min(nR+dx, nR) - dx;
		mCnt((iL:iR) + dx, :) = mCnt((iL:iR) + dx, :) + 1;
		sx((iL:iR) + dx, :) = sx((iL:iR) + dx, :) + img(iL:iR, :);
	end
	sx = sx ./ mCnt;
	mCnt = zeros(nR, nC);
	for dx = -nH:nH,
		iL = max(1+dx, 1) - dx;
		iR = min(nC+dx, nC) - dx;
		mCnt(:, (iL:iR) + dx) = mCnt(:, (iL:iR) + dx) + 1;
		res(:, (iL:iR) + dx) = res(:, (iL:iR) + dx) + sx(:, iL:iR);
	end
	res = res ./ mCnt;
end