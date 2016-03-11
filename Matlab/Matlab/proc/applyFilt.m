function res = applyFilt(imdft, filt, expand)
    imdft = fftshift(imdft, 2);
    siz = size(imdft);
    siz(2) = 1;
    imdft = imdft .* repmat(filt, siz);

    [indices, lenO] = getIDXFromFilter1D(filt, expand);
    res = ifft(ifftshift(imdft(:, indices, :), 2), [], 2) * length(indices);
    res = res / length(filt); 
end
