function depthSelect = getDepthSelect(dsId)
    if dsId == 0 || strcmp(dsId, 'none'),
        depthSelect = 'none';
    elseif dsId == 1 || strcmp(dsId, 'lowD'),
        depthSelect = 'lowD';
    elseif dsId == 2 || strcmp(dsId, 'impD'),
        depthSelect = 'impD';
    elseif dsId == 3 || strcmp(dsId, 'upD'),
        depthSelect = 'upD';
    elseif dsId == 4 || strcmp(dsId, 'gtD'),
        depthSelect = 'gtD';
    elseif dsId == 5 || strcmp(dsId, 'gtImpD'),
        depthSelect = 'gtImpD';
    elseif dsId == 6 || strcmp(dsId, 'impDH'),
        depthSelect = 'impDH';
    elseif dsId == 7 || strcmp(dsId, 'hafD'),
        depthSelect = 'hafD';
    elseif dsId == 8 || strcmp(dsId, 'lowImpD'),
        depthSelect = 'lowImpD';
    elseif dsId == 9 || strcmp(dsId, 'lowFastD')
        depthSelect = 'lowFastD';
    else,
        error('unknown depth select type');
    end
end