function res = appendRef(img, app)
    siz = size(img);
    siz(2) = 1;
    nC = size(img, 2);
    res = circshift(cat(2, img, (img(:, 2*app:-1:1, :) .* repmat((1:2*app)-0.5, siz) +...
            img(:, nC:-1:nC-2*app+1, :) .* repmat((2*app:-1:1)-0.5, siz)) / (2*app)), ...
                    app, 2);
end