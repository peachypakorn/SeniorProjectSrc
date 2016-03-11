% Read image 

slice_v_pos = 49;
%slice_v_pos = 25;

source_image_1 = double(imread('../bbb_106_left.jpg')) / 255;
source_image_2 = double(imread('../bbb_106_right.jpg')) / 255;
%source_image_1 = source_image_1(421:740, 1:1919, :);
%source_image_2 = source_image_2(421:740, 1:1919, :);
source_image_1 = source_image_1(541:620, :, :);
source_image_2 = source_image_2(541:620, :, :);
%source_image_1 = source_image_1(:, 1:1919, :);
%source_image_2 = source_image_2(:, 1:1919, :);

%source_image_1 = double(imread('../bbb_252p7_left.jpg')) / 255;
%source_image_2 = double(imread('../bbb_252p7_right.jpg')) / 255;
%source_image_1 = source_image_1(321:640, 1:1919, :);
%source_image_2 = source_image_2(321:640, 1:1919, :);
%source_image_1 = source_image_1(:, 1:1919, :);
%source_image_2 = source_image_2(:, 1:1919, :);

%source_image_1 = rgb2ntsc(double(imread('../frame1_r.jpg')) / 255 );
%source_image_2 = rgb2ntsc(double(imread('../frame2_r.jpg')) / 255 );

%source_image_1 = double(imread('../transparency_left.jpg')) / 255;
%source_image_2 = double(imread('../transparency_right.jpg')) / 255;

%source_image_1 = double(imread('../left.jpg')) / 255;
%source_image_2 = double(imread('../right.jpg')) / 255;

%source_image_1 = double(imread('../frame1.jpg')) / 255;
%source_image_2 = double(imread('../frame2.jpg')) / 255;
%mask_1 = double(imread('../mask1.jpg')) / 255;
%mask_2 = double(imread('../mask2.jpg')) / 255;
%source_image_1 = source_image_1 .* mask_1;
%source_image_2 = source_image_2 .* mask_2;

%source_image = double(imread('test/synthesis.png')) / 255;
%source_image_1 = source_image(:, 1:size(source_image, 2)/2, :);
%source_image_2 = source_image(:, size(source_image, 2)/2+1:size(source_image, 2), :);

%source = double(imread('../video_1_aligned_16.jpg')) / 255;
%source_image_1 = source(:, 1:size(source, 2)/2, :);
%source_image_2 = source(:, size(source, 2)/2+1:size(source, 2), :);

%source_image_1 = double(imread('../river_001_left.jpg')) / 255;
%source_image_2 = double(imread('../river_001_right.jpg')) / 255;

%% Make Signal
signalWidth = size(source_image_1, 2);
x = linspace(-5,5,signalWidth);
%magnifications = [0 1 2 3];
magnifications = [-0.5];
rho = 0.0;
%% Create filters
filters =  get1DCosFilters(signalWidth, 128);

outimg = zeros([ size(source_image_1)  length(magnifications)*2 ]);
qualGhost = zeros([ size(source_image_1) length(magnifications)*2 ]);
qualDispn = zeros([ size(source_image_1) 2 ]);

[pyr, pind] = buildSCFpyrGen1D(squeeze(source_image_1(1,:,1)), filters, 'alias', true);

pyr = zeros([ size(pyr,1) size(source_image_1, 1)]);
pyr2 = zeros([ size(pyr,1) size(source_image_1, 1)]);

weight = [0.3;0.6;0.1];

for i=1:size(source_image_1, 1),
    %vec1 = squeeze(source_image_1(i,:, :)) * weight;
    %vec2 = squeeze(source_image_2(i,:, :)) * weight;
    vec1 = rgb2ntsc(squeeze(source_image_1(i,:, :)));
    vec2 = rgb2ntsc(squeeze(source_image_2(i,:, :)));
    [pyr(:,i), ~] = buildSCFpyrGen1D(reshape(vec1(:,1), 1, [size(vec1, 1)]), filters, 'alias', true);
    %tmp = buildSCFpyrGen1D(reshape(mask_1(i, :, 1), 1, [size(mask_1, 2)]), filters);
    %tmp(tmp == 0) = 1;
    %pyr(:, i) = pyr(:, i) ./ tmp;
    [pyr2(:,i), ~] = buildSCFpyrGen1D(reshape(vec2(:,1), 1, [size(vec1, 1)]), filters, 'alias', true);
    %tmp = buildSCFpyrGen1D(reshape(mask_2(i, :, 1), 1, [size(mask_1, 2)]), filters);
    %tmp(tmp == 0) = 1;
    %pyr2(:, i) = pyr2(:, i) ./ tmp;
    
end;
%% Computing Filter
nFilts = max(size(filters));
d1 = angle(pyr);
d2 = angle(pyr2);
phase_shift = (d2 - d1);

%% gaussian filter
%std = 10;
%ext = 20;
%extscale = 1;
%h = exp(-(-ext:extscale:ext).^2/(2*std^2)) / (std*sqrt(2*pi));
%for k = 1:nFilts,
%    ind = pyrBandIndices(pind,k);
%    phase_shift( ind, :) = imfilter( exp(phase_shift( ind, :)*1i) , h);
%end


sigma_v = 10;
sigma_p_base = 2.0;
sigma_p = zeros([size(pyr,1) size(source_image_1, 1)]);
for k = 1:nFilts,
    ind = pyrBandIndices(pind,k);
    sigma_p(ind, :) = sigma_p_base * 1.^(1-k);
end
phase_vec = exp(phase_shift*1i);
phase_vec_sum = exp(phase_shift*1i);
phase_vec_weight = ones([ size(pyr,1) size(source_image_1, 1)]);
v = size(source_image_1, 1);
for d = 1:sigma_v*2,
    weight = exp(-(abs(phase_vec(:, 1+d:v) - phase_vec(:, 1:v-d)).^2 ./ (2 * sigma_p(:, 1+d:v).^2)) - ...
                 (d).^2 / (2 * sigma_v^2));
    phase_vec_sum(:, 1:v-d) = phase_vec_sum(:, 1:v-d, :) + phase_vec(:, 1+d:v) .* weight;
    phase_vec_weight(:, 1:v-d) = phase_vec_weight(:, 1:v-d) + weight;
    weight = exp((d).^2 / (2 * sigma_v^2));
    phase_vec_sum(:, v-d+1:v) = phase_vec_sum(:, v-d+1:v, :) + phase_vec(:, v) * weight * ones(1, d);
    phase_vec_weight(:, v-d+1:v) = phase_vec_weight(:, v-d+1:v) + weight;
    weight = exp(-(abs(phase_vec(:, 1:v-d) - phase_vec(:, 1+d:v)).^2 ./ (2 * sigma_p(:, 1+d:v).^2)) - ...
                 (d).^2 / (2 * sigma_v^2));
    phase_vec_sum(:, 1+d:v) = phase_vec_sum(:, 1+d:v, :) + phase_vec(:, 1:v-d) .* weight;
    phase_vec_weight(:, 1+d:v) = phase_vec_weight(:, 1+d:v) + weight;
    weight = exp((d).^2 / (2 * sigma_v^2));
    phase_vec_sum(:, 1:d) = phase_vec_sum(:, 1:d, :) + phase_vec(:, 1) * weight * ones(1, d);
    phase_vec_weight(:, 1:d) = phase_vec_weight(:, 1:d) + weight;
end;
 
phase_shift = angle(phase_vec_sum ./ phase_vec_weight);

multiphasedelta = abs(mod(phase_shift + 5*pi, 2*pi) - pi);

if false,
    slice_vis = zeros([(nFilts+1)*10, size(source_image_1, 2), 3]);
    for k = nFilts:-1:1,
        ind = pyrBandIndices(pind,k);
        for y = 1:size(source_image_1, 2),
            if y > (2^max(0,k-2)) * size(ind, 2),
                slice_vis(k*10+1:k*10+10, y, :) = ones([10, 1]) * [0, 0, 0];
            else,
                phase = 6*multiphasedelta(ind(int32(floor((y-1) / (2^max(0,k-2)))) + 1), slice_v_pos) / (2^(nFilts-2-max(0,k-2)));
                slice_vis(k*10+1:k*10+10, y, :) = ones([10, 1]) * [phase / (2*pi), 0, 0];
            end;
        end;
    end;
    slice_vis = ntsc2rgb(slice_vis);
    slice_vis(1:5, :, :) = source_image_1(slice_v_pos-2:slice_v_pos+2, :, :);
    slice_vis(6:10, :, :) = source_image_2(slice_v_pos-2:slice_v_pos+2, :, :);
end

if false,
    slice_vis2 = zeros([(nFilts+2)*10, size(source_image_1, 2), 3]);
    for k = nFilts:-1:1,
        ind = pyrBandIndices(pind,k);
        for y = 1:size(source_image_1, 2),
            if y > (2^max(0,k-2)) * size(ind, 2),
                slice_vis2(k*10+1:k*10+10, y, :) = ones([10, 1]) * [0.0, 0, 0];
            else,
                phase = 6*multiphasedelta(ind(int32(floor((y-1) / (2^max(0,k-2)))) + 1), slice_v_pos) / (2^(nFilts-2-max(0,k-2)));
                slice_vis2(k*10+1:k*10+10, y, :) = ones([10, 1]) * [phase / (2*pi), 0, 0];
            end;
        end;
    end;
    slice_vis2 = ntsc2rgb(slice_vis2);
    slice_vis2(1:5, :, :) = source_image_1(slice_v_pos-2:slice_v_pos+2, :, :);
    slice_vis2(6:10, :, :) = source_image_2(slice_v_pos-2:slice_v_pos+2, :, :);
end

for i=1:size(source_image_1, 1),
    [pyr(:,i), ~] = buildSCFpyrGen1D(squeeze(source_image_1(i,:, c)), filters, 'alias', true);
    [pyr2(:,i), ~] = buildSCFpyrGen1D(squeeze(source_image_2(i,:, c)), filters, 'alias', true);
end;
%multiphasedelta(pyrBandIndices(pind,1), :) = multiphasedelta(pyrBandIndices(pind,2), :);

%% Reconstructing
for m=1:length(magnifications),
    mag = magnifications(m);

    pyrghost = abs(phase_shift * (mag+1)) * 2 / pi;
    pyrghost = clip(pyrghost - 1, 0, 2);
    pyrghost1 = pyrghost .* (abs(pyr) .^ 2);

    for i=1:size(source_image_1, 1),       
        qualGhost(i, :, :, m*2 - 1) = repmat(reconSCFpyrGen1D(pyrghost1(:,i),pind,filters, 'alias', true, 'smooth', true)', 1, 3);
    end
        
    pyrghost2 = pyrghost .* (abs(pyr2) .^ 2);
    for i=1:size(source_image_1, 1),   
        qualGhost(i, :, :, m*2 ) = repmat(reconSCFpyrGen1D(pyrghost2(:,i),pind,filters, 'alias', true, 'smooth', true)', 1, 3);
    end

end

pyrdispn1 = zeros(size(pyr));
pyrdispn2 = zeros(size(pyr));
for f=1:nFilts,
    start = sum(pind(1:f-1, 2)) + 1;
    len = pind(f, 2);
    if len > 2,
        for j = 0:len-1
            nxt = mod(j+1, len);
            prv = mod(j+len-1, len);
            pyrdispn1(start+j, :) = (((phase_shift(start+j, :) - phase_shift(start+nxt, :)) .^ 2) ...
                                    .* abs(pyr(start+j, :)) .* abs(pyr(start+nxt, :)) + ...
                                    ((phase_shift(start+j, :) - phase_shift(start+prv, :)) .^ 2) ...
                                    .* abs(pyr(start+j, :)) .* abs(pyr(start+nxt, :))) / (pi^2);
            pyrdispn2(start+j, :) = (((phase_shift(start+j, :) - phase_shift(start+nxt, :)) .^ 2) ...
                                    .* abs(pyr2(start+j, :)) .* abs(pyr2(start+nxt, :)) + ...
                                    ((phase_shift(start+j, :) - phase_shift(start+prv, :)) .^ 2) ...
                                    .* abs(pyr2(start+j, :)) .* abs(pyr2(start+nxt, :))) / (pi^2);
        end
    end
end

for i=1:size(source_image_1, 1),
    qualDispn(i, :, :, 1) = repmat(reconSCFpyrGen1D(pyrdispn1(:,i),pind,filters, 'alias', true, 'smooth', true)', 1, 3);
    qualDispn(i, :, :, 2) = repmat(reconSCFpyrGen1D(pyrdispn2(:,i),pind,filters, 'alias', true, 'smooth', true)', 1, 3);
end

for c=1:size(source_image_1,3),
    for i=1:size(source_image_1, 1),
        [pyr(:,i), ~] = buildSCFpyrGen1D(squeeze(source_image_1(i,:, c)), filters, 'alias', true);
        [pyr2(:,i), ~] = buildSCFpyrGen1D(squeeze(source_image_2(i,:, c)), filters, 'alias', true);
    end;
    %multiphasedelta(pyrBandIndices(pind,1), :) = multiphasedelta(pyrBandIndices(pind,2), :);
    
    %% Reconstructing
    for m=1:length(magnifications),
        mag = magnifications(m);

        pyrnew = exp( -( multiphasedelta * rho ).^2 / 2).* exp( 1i* (-phase_shift) * mag) .* pyr;
        for i=1:size(source_image_1, 1),       
            recons = reconSCFpyrGen1D(pyrnew(:,i),pind,filters, 'alias', true) ;
            outimg(i, :, c, m*2 - 1) = recons;
        end
            
        pyrnew = exp( -( multiphasedelta * rho ).^2 / 2).* exp( 1i* (phase_shift * (mag))) .* pyr2;
        for i=1:size(source_image_1, 1),       
            recons = reconSCFpyrGen1D(pyrnew(:,i),pind,filters, 'alias', true) ;
            outimg(i, :, c, m*2 ) = recons;
        end
%        plot(x, signal,'r') ; hold on;plot(x, signal2,'g') ; plot(x, recons,'b') ; hold off

    end
end

if false,
    slice_vis2(nFilts*10+11:nFilts*10+15, :, :) = repmat(outimg(slice_v_pos,:,:, 7), [5, 1, 1]);
    slice_vis2(nFilts*10+16:nFilts*10+20, :, :) = repmat(outimg(slice_v_pos,:,:, 8), [5, 1, 1]);
end

out_folder = sprintf('bbb_106_midpoint_%.2f', rho);
mkdir(out_folder);
for m=1:length(magnifications),
    outimg(:,:,:, 2*m-1) = outimg(:,:,:, 2*m-1);
    outimg(:,:,:, 2*m) = outimg(:,:,:, 2*m);
    imwrite(outimg(:,:,:, 2*m-1), sprintf('%s/%03d.png', out_folder, length(magnifications)-m+1),'png');
    imwrite(outimg(:,:,:, 2*m), sprintf('%s/%03d.png', out_folder, length(magnifications)+m),'png');
    imwrite(qualGhost(:,:,:, 2*m-1)*10, sprintf('%s/qg-%03d.png', out_folder, length(magnifications)-m+1),'png');
    imwrite(qualGhost(:,:,:, 2*m)*10, sprintf('%s/qg-%03d.png', out_folder, length(magnifications)+m),'png');
end
imwrite(qualDispn(:,:,:, 1)*10, sprintf('%s/qd-1.png', out_folder),'png');
imwrite(qualDispn(:,:,:, 2)*10, sprintf('%s/qd-2.png', out_folder),'png');
%imwrite(slice_vis, sprintf('%s/vis_1.png', out_folder),'png');
%imwrite(slice_vis2, sprintf('%s/vis_2.png', out_folder),'png');

