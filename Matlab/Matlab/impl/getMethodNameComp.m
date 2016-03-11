function name = getMethodNameComp(methodID)
    if methodID == 1 || strcmp(methodID, 'eol'),
        name = 'eol';
    elseif methodID == 2 || strcmp(methodID, 'lagr'),
        name = 'lagr';
    elseif methodID == 3 || strcmp(methodID, 'eol2'),
        name = 'eol2';
    elseif methodID == 4 || strcmp(methodID, 'lagrWav'),
        name = 'lagrWav';
    elseif methodID == 5 || strcmp(methodID, 'eulaProp'),
        name = 'eulaProp';
    elseif methodID == 6 || strcmp(methodID, 'eulaPropHf'),
        name = 'eulaPropHf';
    elseif methodID == 7 || strcmp(methodID, 'lagrSmooth'),
        name = 'lagrSmooth';
    else,
        error('unknown method');
    end
end