function [outimgL, outimgR] = viewSimulate(imgs, view_pos, view_dist)
    pixel_size = 0.1989;
    view_angle = 20;
    eye_sep_original = 65;
    opt_view_dist = 400;

    width = size(imgs, 2);
    height= size(imgs, 1);
    nCol = size(imgs, 3);
    num_views  = size(imgs, 4);

    subpixel_size = pixel_size/num_views;
    eye_sep = opt_view_dist*tan(view_angle/2/180*pi)*2/num_views;
    opt_spacer = subpixel_size * opt_view_dist / eye_sep;
    mask_cycle = subpixel_size/(1+subpixel_size/eye_sep)*num_views;

    for vpos = view_pos + width/2 + [-eye_sep_original/2, eye_sep_original/2],
        i = repmat(1:width, [height, 1]);
        j = repmat((1:height)', [1, width]);
        x = (i-j+height) * mask_cycle;
        pixel_x = x-(vpos-x)*opt_spacer/view_dist;
        pixel_x = pixel_x - floor(pixel_x/pixel_size)*pixel_size+pixel_size;
        pixel_i = floor(pixel_x/subpixel_size);
        w = (pixel_x-subpixel_size*pixel_i)/subpixel_size;
        pixel_i = mod(pixel_i+num_views, num_views) + 1;
        pixel_i2 = mod(pixel_i-2+num_views, num_views) + 1;

        q00 = imgs(sub2ind([height, width, nCol, num_views], repmat(j, [1, 1, nCol]), ...
                   repmat(i, [1, 1, nCol]), repmat(reshape(1:nCol, [1, 1, nCol]), [height, width, 1]), ...
                   repmat(pixel_i, [1, 1, nCol])));
        q01 = imgs(sub2ind([height, width, nCol, num_views], repmat(j, [1, 1, nCol]), ...
                   repmat(i, [1, 1, nCol]), repmat(reshape(1:nCol, [1, 1, nCol]), [height, width, 1]), ...
                   repmat(pixel_i2, [1, 1, nCol])));
        result = repmat(w, [1, 1, nCol]) .* q00 + repmat((1-w), [1, 1, nCol]) .* q01;
        if vpos < view_pos + width/2,
            outimgL = result;
        else,
            outimgR = result;
        end
    end
end