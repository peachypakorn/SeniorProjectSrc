w = 800;
h = 400;
back_cell = 20;
front_cell = 5;
st = [256 57 ];
nc = [60 60];
ex = 15;
back_col = [0.7 0.7 0.7];
front_col = [0 0 0];

output = zeros([h w 3 2]);


for ii = 1:2,
    figure;
    imshow(output(:,:,:,ii));
    hold on;
    set(gca,'Units','normalized','Position',[0 0 1 1]);  %# Modify axes size
    set(gcf,'Units','pixels','Position',[200 200 w h]);  %# Modify figure size
  
    for wi = 1:back_cell:w,
        for hi = 1:back_cell:h,
            if mod(wi+hi-2, 2*back_cell) == back_cell,
                rectangle('Position', [wi hi back_cell back_cell], 'FaceColor', [1 1 1], 'LineStyle', 'none');
            else,
                rectangle('Position', [wi hi back_cell back_cell], 'FaceColor', back_col, 'LineStyle', 'none');
            end
        end
    end
    f = getframe(gcf);              %# Capture the current window
    output(:,:,:,ii)= f.cdata;
    close;
    figure;
    imshow(output(:,:,:,ii)/256);
    hold on;
    set(gca,'Units','normalized','Position',[0 0 1 1]);  %# Modify axes size
    set(gcf,'Units','pixels','Position',[200 200 w h]);  %# Modify figure size

    for wi = 1:nc(1),
        for hi = 1:nc(2),
            p1 = st + [wi-1 hi-1 ] .* [front_cell+2*(ii-1)*ex*(hi-1)/(nc(1)*nc(2)) front_cell ] - [ex*(ii-1)*(hi-1)/nc(2) 0 ];
            p2 = st + [wi hi-1 ] .* [front_cell+2*(ii-1)*ex*(hi-1)/(nc(1)*nc(2)) front_cell ] - [ex*(ii-1)*(hi-1)/nc(2) 0 ];
            p3 = st + [wi hi ] .* [front_cell+2*(ii-1)*ex*(hi)/(nc(1)*nc(2)) front_cell ] - [ex*(ii-1)*(hi)/nc(2) 0 ];
            p4 = st + [wi-1 hi ] .* [front_cell+2*(ii-1)*ex*(hi)/(nc(1)*nc(2)) front_cell ] - [ex*(ii-1)*(hi)/nc(2) 0 ];
            P = [p1;p2;p3;p4];
            if mod(wi+hi, 2) == 1,
                fill(P(:,1), P(:,2), [1 1 1],'FaceColor', [1 1 1], 'LineStyle', 'none');
            else,
                fill(P(:,1), P(:,2), front_col, 'FaceColor', front_col, 'LineStyle', 'none');
            end
        end
    end
    myaa(4);
    f = getframe(gcf);              %# Capture the current window
    output(:,:,:,ii)= f.cdata;  %# Save the frame data
    close;
    close;
end


out_folder = sprintf('data/');

imwrite(output(:, :, :, 1)/255, sprintf('%s/check9_1.png', out_folder),'png');
imwrite(output(:, :, :, 2)/255, sprintf('%s/check9_2.png', out_folder),'png');
