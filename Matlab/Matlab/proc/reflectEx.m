function reflect = reflectEx(img, nEx)
    if isscalar(img)
        reflect = img;
    else
        nR = size(img, 1);
        nC = size(img, 2);
        if nargin < 2,
            nEx = round(nC * 0.1);
        end
        nCh = size(img, 3);
        reflect = zeros(nR, nC+2*nEx, nCh);
        reflect(:, nEx+(1:nC), :) = img;
        reflect(:, 1:nEx, :) = img(:, nEx:-1:1, :);
        reflect(:, nC+nEx+(1:nEx), :) = img(:, nC+1-(1:nEx), :);
    end
end