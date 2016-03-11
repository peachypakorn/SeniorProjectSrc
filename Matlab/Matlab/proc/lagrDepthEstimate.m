
function [depthL, depthR] = lagrDepthEstimate(Il, Ir, maxDLagr, dReg)
    % https://www.ims.tuwien.ac.at/publications/tuw-210567

    % Parameter settings
    r = 9;                  % filter kernel in eq. (3) has size r \times r
    eps = 0.0001;           % \epsilon in eq. (3)
    thresColor = 7/255;     % \tau_1 in eq. (5)
    thresGrad = 2/255;      % \tau_2 in eq. (5)
    gamma = 0.11;           % (1- \alpha) in eq. (5)
    threshBorder = 3/255;   % some threshold for border pixels
    gamma_c = 0.1;          % \sigma_c in eq. (6)
    gamma_d = 9;            % \sigma_s in eq. (6)
    r_median = 19;          % filter kernel of weighted median in eq. (6) has size r_median \times r_median

    % Convert to grayscale
    Il_g = rgb2gray(Il);
    Ir_g = rgb2gray(Ir);

    % Mirror image
    Il_1 = flipdim(Ir,2);
    Ir_1 = flipdim(Il,2);

    % Compute gradient in X-direction from grayscale images
    fx_l = gradient(Il_g);
    fx_r = gradient(Ir_g);
    fx_l = fx_l+0.5; % To get a range of values between 0 to 1
    fx_r = fx_r+0.5; % To get a range of values between 0 to 1
    fx_l_1 = flipdim(fx_r,2);
    fx_r_1 = flipdim(fx_l,2);

    [m,n,c] = size(Il);

    dispVol = ones(m,n,2*maxDLagr+1)*(threshBorder);
    dispVol1 = ones(m,n,2*maxDLagr+1)*(threshBorder);

    % Create initial cost volume (eq. (5) in the paper)
    for d=-maxDLagr:maxDLagr
        if d > 0,
            range1 = d+1:n;
            range2 = 1:n-d;
        else,
            range1 = 1:n+d;
            range2 = -d+1:n;
        end

        % Truncated SAD of color images for current displacement
        tmp = ones(m,n,c)*threshBorder;
        tmp(:,range1,:) = Ir(:,range2,:);
        p_color = abs(tmp - Il);
        p_color = sum(p_color,3)*0.333333333333;
        p_color = min(p_color,thresColor);

        % Truncated SAD of gradient images for current displacement
        tmp = ones(m,n)*threshBorder;
        tmp(:,range1) = fx_r(:,range2);
        p_grad = abs(tmp - fx_l);
        p_grad = min(p_grad,thresGrad);

        p = gamma*p_color + (1-gamma)*p_grad; % Combined color and gradient

        % Same for other view
        tmp1 = ones(m,n,c)*threshBorder;
        tmp1(:,range1,:) = Ir_1(:,range2,:);
        p1_color = abs(tmp1 - Il_1);
        p1_color = sum(p1_color,3)*0.333333333333;
        p1_color = min(p1_color,thresColor);

        tmp1 = ones(m,n)*threshBorder;
        tmp1(:,range1) = fx_r_1(:,range2);
        p1_grad = abs(tmp1 - fx_l_1);
        p1_grad = min(p1_grad,thresGrad);

        p1 = gamma*p1_color + (1-gamma)*p1_grad; % Combined color and gradient

        % Set value in cost volume
        dispVol(:,:,d+1+maxDLagr) = p;
        dispVol1(:,:,d+1+maxDLagr) = flipdim(p1,2);
        %imshow(dispVol(:,:,d+1+maxDLagr) * 50);
        %figure;
    end

    % Smooth cost volume with guided filter (using eq. (1) with weights in (4))
    for d=1:2*maxDLagr+1 % use regular for loop when not using the parallel computing toolbox
        p = dispVol(:,:,d);
        p1 = dispVol1(:,:,d);

        q = guidedfilter_color(Il, double(p), r, eps);
        p1 =  flipdim(p1,2);
        q1 = guidedfilter_color(Il_1, double(p1), r, eps);

        dispVol(:,:,d) = q * (1 + ((d - maxDLagr - 1) / dReg) ^2);
        dispVol1(:,:,d) = flipdim(q1,2) * (1 + ((d - maxDLagr - 1) / dReg) ^2);
    end

    % Winner take all label selection (eq. (2))
    [unused,labels_left] = min(dispVol,[],3);
    [unused,labels_right] = min(dispVol1,[],3);
    labels_left = labels_left - (1+maxDLagr);
    labels_right = labels_right - (1+maxDLagr);

    %figure
    %ana = zeros([size(labels_left), 3]);
    %ana(:, :, 1) = 0.5*labels_left / maxDLagr + 0.5;
    %ana(:, :, 2:3) = repmat(0.5*labels_right / maxDLagr + 0.5, [1, 1, 2]);
    %imshow(ana);

    % Left-right consistency check
    Y = repmat((1:m)', [1 n]);
    X = repmat(1:n, [m 1]) - labels_left;
    X(X<1) = 1;
    X(X>n) = n;
    indices = sub2ind([m,n],Y,X);

    %imshow(labels_left / (2*maxDLagr) + 0.5);
    %figure;
    %imshow(labels_right / (2*maxDLagr) + 0.5);
    %figure;
    %pause
    final_labels_left = labels_left;
    final_labels_left(abs(labels_left - labels_right(indices))>=1) = -maxDLagr-1;

    %figure
    %ana(:, :, 1) = 0.5*final_labels_left / maxDLagr + 0.5;

    % Fill and filter (post-process) pixels that fail the consistency check
    final_labels_left = fillPixelsReference(Il, final_labels_left, gamma_c, gamma_d, r_median, maxDLagr);

    depthL = final_labels_left;

    % Left-right consistency check
    Y = repmat((1:m)', [1 n]);
    X = repmat(1:n, [m 1]) + labels_right;
    X(X<1) = 1;
    X(X>n) = n;
    indices = sub2ind([m,n],Y,X);

    final_labels_right = labels_right;
    final_labels_right(abs(labels_right - labels_left(indices))>=1) = -maxDLagr-1;

    %ana(:, :, 2:3) = repmat(0.5*final_labels_right / maxDLagr + 0.5, [1, 1, 2]);
    %imshow(ana);

    % Fill and filter (post-process) pixels that fail the consistency check
    final_labels_right = fillPixelsReference(Ir, final_labels_right, gamma_c, gamma_d, r_median, maxDLagr);

    depthR = final_labels_right;
end
