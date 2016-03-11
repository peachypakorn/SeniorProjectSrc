function name = getMethodName(methodID)
    if methodID == 1 || strcmp(methodID, 'eulaOnLagr'),
        name = 'eulaOnLagr';
    elseif methodID == 2 || strcmp(methodID, 'lagrangian'),
        name = 'lagrangian';
    elseif methodID == 3 || strcmp(methodID, 'eularian2D'),
        name = 'eularian2D';
    elseif methodID == 4 || strcmp(methodID, 'eulaSynth'),
        name = 'eulaSynth';
    elseif methodID == 5 || strcmp(methodID, 'eulaOnLagrRefl'),
        name = 'eulaOnLagrRefl';
    else,
        error('unknown method');
    end
end