function res = eularian2D(img1, img2, testParam, runParam)
    rho = runParam.rho2D;
    ht = runParam.ht;
    nA = runParam.nA;
    dA = runParam.dA;
    if isfield(runParam, 'remap'),
        maxDepth = runParam.remapDepth;
    else
        maxDepth = 0;
    end
    if isfield(runParam, 'orients'),
        orients = runParam.orients;
    else
        orients = 8;
    end
    twidth = 1;

    nLevels  = maxSCFpyrHt(img1(:, :, 1));
    img1 = rgb2ntsc(img1);
    img2 = rgb2ntsc(img2);

    [h w nC] = size(img1);

    res = zeros(h,w,nC,2*nA);

    filt = orients-1;
    [~, pindO] = buildSCFpyr(img1(:,:,1), nLevels, filt);
    numScales = (size(pindO,1)-2)/orients + 2;
    numBands = size(pindO,1);
    numElements = dot(pindO(:,1),pindO(:,2));

    %% Scale up magnification levels
    %  We use alpha 1 and manipulate directly phase
    magPhase = [0 repmat(1, [1, nLevels]) 0]';
    magPhase = scaleBand2all(magPhase, numScales, orients);

    %% Create complex steerable pyramids
    pyrCA1 = buildSCFpyrC(img1, nLevels, filt, twidth);
    pyrCA2 = buildSCFpyrC(img2, nLevels, filt, twidth);
	disp('Pyramid Bulit');

    % fixing phase shifts to be accurate beyond 2pi range

    sizes = pindO(:,1).*pindO(:,2);
    c = [0; sizes];
    c = cumsum(c);

    phase_shift = angle(pyrCA2) - angle(pyrCA1);

    for level = (size(pindO,1) - 9):-1:2
        previous_phases = phase_shift(c(level+8)+1:c(level+8+1));
        previous_phases = reshape(previous_phases, pindO(level+8,:));
        previous_phases = imresize(previous_phases, pindO(level,:), 'nearest');

        current_phases = phase_shift(c(level)+1:c(level+1));
        current_phases = reshape(current_phases, pindO(level,:));

        %mask = pi < abs(2*previous_phases);
        %current_phases = max(current_phases, 2* previous_phases);
        %current_phases(mask) = previous_phases(mask);
        %phase_shift(c(level)+1:c(level+1)) = current_phases;
        %imshow(current_phases / (2*pi) + 0.5);
        %figure;
        current_phases = current_phases + round((2 * previous_phases - current_phases) / (2*pi)) * 2 * pi;
        phase_shift(c(level)+1:c(level+1)) = current_phases;
    end

    aa_phase = abs(phase_shift);
    phase_shift_d = phase_shift;

    for level = (size(pindO,1) - 9):-1:2
        previous_ps = phase_shift_d(c(level+8)+1:c(level+8+1));
        previous_ps = reshape(previous_ps, pindO(level+8,:));
        previous_ps = imresize(previous_ps, pindO(level,:), 'nearest');
        current_ps = phase_shift_d(c(level)+1:c(level+1));
        current_ps = reshape(current_ps, pindO(level,:));

        previous_phases = aa_phase(c(level+8)+1:c(level+8+1));
        previous_phases = reshape(previous_phases, pindO(level+8,:));
        previous_phases = imresize(previous_phases, pindO(level,:), 'nearest');

        current_phases = aa_phase(c(level)+1:c(level+1));
        current_phases = reshape(current_phases, pindO(level,:));

        current_ps(current_phases < 2*previous_phases) = 2*previous_ps(current_phases < 2*previous_phases);
        current_phases = max(current_phases, 2*previous_phases);
        current_phases = current_phases(:);
        aa_phase(c(level)+1:c(level+1)) = current_phases;
        phase_shift_d(c(level)+1:c(level+1)) = reshape(current_ps, [], 1);
    end

    for i=2:9
        previous_ps = phase_shift_d(c(i)+1:c(i+1));
        previous_ps = reshape(previous_ps, pindO(i,:));
        previous_ps = imresize(previous_ps, pindO(1,:), 'nearest');
        current_ps = phase_shift_d(c(1)+1:c(1+1));

        previous_phases = aa_phase(c(i)+1:c(i+1));
        previous_phases = reshape(previous_phases, pindO(i,:));
        previous_phases = imresize(previous_phases, pindO(1,:), 'nearest');

        current_phases = aa_phase(c(1)+1:c(1+1));

        current_ps(current_phases < previous_phases(:)') = previous_ps(current_phases < previous_phases(:)');
        aa_phase(c(1)+1:c(1+1)) = max(current_phases, previous_phases(:)');
        phase_shift_d(c(1)+1:c(1+1)) = current_ps;
    end

    if maxDepth > 0,
        fact = 4*ones(size([[phase_shift_d]])) / (2*pi);
        for i = 2:size(pindO,1) - 1
            i1 = floor((i-2) / 8);
            i2 = i-2-i1*8;
            fact(c(i)+1:c(i+1)) = fact(c(i)+1:c(i+1)) * (2^i1)*2 / (0.5 + (sin(2*(i2+1)*pi / 8) - sin(2*i2*pi / 8)) / (pi/2));
        end
        phase_clip = maxDepth * atan(phase_shift .* fact / maxDepth) ./ fact;
        aa_clip = maxDepth * atan(aa_phase .* fact / maxDepth) ./ fact;
    else
        phase_clip = phase_shift;
        aa_clip = aa_phase;
    end
	disp('Phase Processed');

    for iA = 1:nA
        phase_mov = (iA - 0.5) * dA * phase_clip - 0.5 * phase_shift;
        res(:,:,:,nA+1-iA) = ntsc2rgb(reconSCFpyrC(exp( -( dA * abs(aa_clip) * rho).^2 / 2) .* ...
                exp(1i*(-phase_mov)) .* pyrCA1, pindO, 'all', 'all', twidth));
        res(:,:,:,nA+iA) = ntsc2rgb(reconSCFpyrC(exp( -( dA * abs(aa_clip) * rho).^2 / 2) .* ...
                exp(1i*(phase_mov)) .* pyrCA2, pindO, 'all', 'all', twidth));
    end
end
