%gtpaths = {'~/Downloads/Aloe/view2.png', '~/Downloads/Dolls/view2.png', '~/Downloads/Art/view2.png'};
%testnames = {'aloe', 'doll', 'art'};
%cuts = [-30, -30, -30];

%gtpaths = {'~/Downloads/Reindeer/view2.png', '~/Downloads/Aloe/view2.png'};
%testnames = {'reindeer', 'aloe'};
%cuts = [-30, -30];

gtpaths = {'./data/testBox6-hd-R2.png'};
testnames = {'testBox6'};
cuts = [0];

%metnames = {'eol-lowD', 'lagr-impD', 'eol-gtD'};
%metnames = {'lagr-gtD', 'lagr-gtImpD'};
metnames = {'eol-lowD', 'lagr-impD', 'eol-gtD', 'lagr-gtD', 'lagr-gtImpD'};

for tid = 1:length(gtpaths),
	for mid = 1:length(metnames),
		img1 = double(imread(gtpaths{tid})) / 255;
		folder = sprintf('./result/compare/%s/plain/%s', testnames{tid}, metnames{mid});
		img2 = double(imread(sprintf('%s/001_003.jpg', folder))) / 255;

		cut = cuts(tid);
		nC = size(img1, 2);
		nR = size(img1, 1);
		if cut < 0,
		    img1 = img1(:, 1:nC+2*cut, :);
		    img2 = imresize(img2, [nR, nC+cut]);
		    img2 = img2(:, 1-cut:nC+cut, :);
		elseif cut > 0,
		    img1 = img1(:, 1+2*cut:nC, :);
		    img2 = imresize(img2, [nR, nC-cut]);
		    img2 = img2(:, 1:nC-2*cut, :);
		end

		%mwrite(0.5 * (1 - structureSimilarity(img1, img2)), sprintf('%s/dssim_003.jpg', folder));
		imwrite(0.5 * (1 - structureSimilarityMulti(img1, img2)), sprintf('%s/dssimM_003.jpg', folder));
	end
end