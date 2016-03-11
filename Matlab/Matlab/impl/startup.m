!synclient HorizTwoFingerScroll=0
show = @(im, kf) {figure(kf), imshow(im, 'Border', 'tight')};
resize = @(im, ref) imresize(im, [size(ref, 1), size(ref, 2)], 'nearest');
