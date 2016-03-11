function simpleFlowSynth1
    amp = [0 1 2 3];
    nA = length(amp);
    trim = 4;
    utrim = 32;
    rho = 1.00;
    rhor = 0.0;
    source_image_1 = evalin('base','exp_image_1');
    depthPyr = evalin('base','depthPyr');
    outimg = zeros([ size(source_image_1) 2*nA ]);
    synth = zeros([ size(source_image_1) nA ]);
    %debug only
    rStripDebug = 16;
    depthFactor = 32;
    depthOffset = 0.5;

    depthImgPyr = cell(length(depthPyr), 1);
    for iLev = 1:length(depthPyr),
        depthImgPyr{length(depthPyr) + 1 - iLev} = imresize(depthPyr{iLev}, [size(source_image_1, 1), size(depthPyr{iLev}, 2)]) * 2^(iLev-1);
    end
    'pyrimid'
    t0 = cputime();
    for r = 1:size(source_image_1, 1),
        for c = 1:3,
            synth(r, :, c, :) = mergableMultiPyrFast(source_image_1(r, :, c), depthImgPyr, r, trim, rho, rhor, amp, utrim);
            %[synth(r, :, c, :), dep] = mergableMultiPyrRec(img(r, :, c), movPyr, r, trim, rho,...
            %        amp, 1, 1, 1, size(img, 2), ones([1 size(img, 2)]), 0);
            %for l = 1:length(dep),
            %    if l < length(dep),
            %       depthP(r, :, :, length(dep)+1-l) = repmat(dep{l}' - dep{l+1}', 1, 3);
            %    else,
            %       depthP(r, :, :, length(dep)+1-l) = repmat(dep{l}', 1, 3);
            %    end
            %end
            %depth(r, :, :) = repmat(dep{1}', 1, 3);
        end
        %r
        if rStripDebug > 0,
            if mod(r, rStripDebug) == 0,
                r
                %for l = 1:length(dep),
                %    imwrite(depthP(:, :, :, l) / min(dscale, size(img, 1) / 2^(l-1)) + 0.5, sprintf('%s/depth%03d.png', out_folder, l),'png');
                %end
                %imwrite(depth(:, :, :) / dscale + doffset, sprintf('%s/depth.png', out_folder),'png');
                for m=1:length(amp),
                    imwrite(repmat(synth(:, :, :, m), [1, 1, 1]), sprintf('%s/%03d.png', 'result/tmp_debug', m+2),'png');
                end
            end
        end
    end
    cputime() - t0

    out_folder = sprintf('result/syn_spaDep_%.2f', rho);
    if rhor > 0,
        out_folder = sprintf('%s_%.2f', out_folder, rhor);
    end

    mkdir(out_folder);
    for m=1:nA,
        imwrite(synth(:, :, :, m), sprintf('%s/%03d.png', out_folder, m),'png');
    end
    for iLev = 1:length(depthPyr),
        imwrite(depthPyr{iLev}*(2^(iLev - 1)) / depthFactor + depthOffset, sprintf('%s/depth-%d.png', out_folder, iLev),'png');
    end
end

function synth = mergableMultiPyrFast(sig, movPyr, r, trim, rho, rhor, amp, utrim)
    nC = length(sig);
    nA = length(amp);
    nLay = round(log2(utrim / trim)) + 1;
    synth = zeros(nC, nA);
    for cLay = 1:nLay,
        cLen = utrim * 2 / (2^cLay);
        sPyr{cLay} = zeros((nC/utrim)*(2^cLay) - 1, cLen, nA);
    end

    for cLay = nLay:-1:1,
        cLen = utrim * 2 / (2^cLay);
        filt = 1:cLen;
        filt = cos((filt - (cLen+1)/2) * pi / cLen);
        lFilt = 1:cLen*2;
        lFilt = cos((lFilt - (2*cLen+1)/2) * pi / (2*cLen));

        dif = 0:cLen-1;
        dif(cLen/2+1) = 0;
        dif(cLen:-1:cLen/2+2) = -1:-1:-cLen/2+1;

        if cLay ~= nLay && cLay ~= 1 % only an estimation
            fac = abs(dif);
            fac(4:2:cLen-2) = 2*fac(4:2:cLen-2);
            fac2 = dif .^ 2 - fac;
        elseif cLay == nLay,
            fac = abs(dif);
            fac(4:2:cLen-2) = 2*fac(4:2:cLen-2);
            fac2 = 0;
        elseif cLay == 1,
            fac = abs(dif);
            fac(4:2:cLen-2) = 2 * fac(4:2:cLen-2);
            fac2 = 0;
        end

        for c = 1:(nC/utrim)*(2^cLay)-1,
            cSt = (c - 1) * cLen / 2 + 1;
            iFilt1 = 1;
            iFilt2 = 1;
            %if cLay == 1,
            %    cFilt = ones(size(filt));
            %else
                if c == 1,
                cFilt = ones(size(filt));
                cFilt(cLen/2+1:cLen) = filt(cLen/2+1:cLen);
            elseif c == (nC/utrim)*(2^cLay)-1,
                cFilt = ones(size(filt));
                cFilt(1:cLen/2) = filt(1:cLen/2);
            else,
                cFilt = filt;
                iFilt1 = ones(size(filt));
                iFilt1(1:cLen/2) = filt(1:cLen/2);
                iFilt2 = ones(size(filt));
                iFilt2(cLen/2+1:cLen) = filt(cLen/2+1:cLen);
            end

            if cLay == nLay,
                cSt = (c - 1) * cLen / 2 + 1;
                sPyr{cLay}(c, :, :) = repmat(sig(cSt:cSt+cLen-1) .* cFilt, [1, 1, nA]);
            end
            pars = [floor(c / 2 + 1e-10), floor((c + 1) / 2 + 1e-10)];
            if pars(1) == pars(2) || cLay == 1,
                pars = [pars(1)];
            end

            for pId = 1:length(pars),
                p = pars(pId);
                if cLay == 1,
                    mov = 0;
                elseif p <= 0 || p >(nC/utrim)*(2^(cLay - 1)) - 1,
                    continue;
                else,
                    mov = (movPyr{cLay - 1}(r, p) + movPyr{cLay - 1}(r, p+1))/2;
                end
                %if cLay == 2,
                %    pFilt = ones(size(lFilt));
                %else
                if p == 1,
                    pFilt = ones(size(lFilt));
                    pFilt(cLen+1:2*cLen) = lFilt(cLen+1:2*cLen);
                elseif p == (nC/utrim)*(2^(cLay-1))-1,
                    pFilt = ones(size(lFilt));
                    pFilt(1:cLen) = lFilt(1:cLen);
                else,
                    pFilt = lFilt;
                end
                cMov = (movPyr{cLay}(r, c) + movPyr{cLay}(r, c+1))/2;
                nMov = cMov - mov;
                for id = 1:nA,
                    if cLay == 1,
                        cSig = sPyr{cLay}(c, :, id);
                        aFilt = cFilt;
                    elseif p * 2 - 1 == c,
                        cSig = sPyr{cLay}(c, :, id) .* pFilt(1:cLen) ./ iFilt1;
                        aFilt = cFilt .* pFilt(1:cLen) ./ iFilt1;
                    elseif p * 2 == c,
                        cSig = sPyr{cLay}(c, :, id) .* pFilt(cLen/2+1:3*cLen/2);
                        aFilt = cFilt .* pFilt(cLen/2+1:3*cLen/2);
                    else,
                        cSig = sPyr{cLay}(c, :, id) .* pFilt(cLen+1:2*cLen) ./ iFilt2;
                        aFilt = cFilt .* pFilt(cLen+1:2*cLen) ./ iFilt2;
                    end
                    avg = sum(cSig) / sum(aFilt);
                    shifted = cSig - avg * aFilt;
                    antialias = exp( -(( 2*pi /cLen )^2) * ((cMov^2)*fac*rho^2 + (nMov^2)*fac2*amp(id)^2*rhor^2) / 2);
                    %if cLay == nLay,
                    %    mMov = movPyr{cLay}(r, c)^2;
                    %    if c > 1,
                    %        mMov = max(mMov, movPyr{cLay}(r, c-1)^2);
                    %    end
                    %    if c < 2^cLay-1,
                    %        mMov = max(mMov, movPyr{cLay}(r, c+1)^2);
                    %    end
                    %    antialias = antialias .* exp( -(( 2*pi /cLen )^2) * (mMov*(dif .^2 - fac)) / 2);
                    %end
                    frq = fft(shifted) .* exp(2*pi*1i*amp(id)*nMov*dif/cLen) .* antialias;
                    shifted = real(ifft(frq));
                    shifted = shifted + avg * aFilt;

                    if cLay ~= 1,
                        if p * 2 - 1 == c,
                            sPyr{cLay - 1}(p, 1:cLen, id) = ...
                                sPyr{cLay - 1}(p, 1:cLen, id) + shifted .* cFilt ./ iFilt1;
                        elseif p * 2 == c,
                            sPyr{cLay - 1}(p, cLen/2+1:3*cLen/2, id) = ...
                                sPyr{cLay - 1}(p, cLen/2+1:3*cLen/2, id) + shifted .* cFilt;
                        else,
                            sPyr{cLay - 1}(p, cLen+1:2*cLen, id) = ...
                                sPyr{cLay - 1}(p, cLen+1:2*cLen, id) + shifted .* cFilt ./ iFilt2;
                        end
                    else,
                        %synth(:, id) = (shifted ./ cFilt)';
                        cSt = (c - 1) * cLen / 2 + 1;
                        synth(cSt:cSt+cLen-1, id) = synth(cSt:cSt+cLen-1, id) + (shifted .* cFilt)';
                    end
                end
            end
        end
    end
end
