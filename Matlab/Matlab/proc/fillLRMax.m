function res = fillLRMin(dep, zero)
	nR = size(dep, 1);
	nC = size(dep, 2);
	res = zeros(nR, nC);
	line = zero * ones(nR, 1);
	for iC = 1:nC,
		mat = dep(:, iC);
		mask = mat ~= zero;
		line(mask) = mat(mask);
		res(:, iC) = line;
	end
	line = zero * ones(nR, 1);
	for iC = nC:-1:1,
		mat = dep(:, iC);
		mask = mat ~= zero;
		line(mask) = mat(mask);
		res(:, iC) = max(res(:, iC), line);
	end
end