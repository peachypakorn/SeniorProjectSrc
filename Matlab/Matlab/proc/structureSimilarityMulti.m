function res = structureSimilarityMulti(img1, img2)
	s = size(img1);
	s = s(1:2);
	res = zeros(s);
	nLayers = 6;
	for iL = 1:nLayers,
		res = res + imresize(structureSimilarity(img1, img2), s) / nLayers;
		img1 = imresize(img1, 0.5);
		img2 = imresize(img2, 0.5);
	end
end