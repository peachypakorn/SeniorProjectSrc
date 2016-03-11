function [q, ax, bx] = guidedfilter_color(I, p, r, eps, weight)
%   GUIDEDFILTER_COLOR   O(1) time implementation of guided filter using a color image as the guidance.
%
%   - guidance image: I (should be a color (RGB) image)
%   - filtering input image: p (should be a gray-scale/single channel image)
%   - local window radius: r
%   - regularization parameter: eps
if nargin <= 4,
    weight = 1;
end

[hei, wid, nP] = size(p);
N = boxfilter(ones(hei, wid), r); % the size of each local patch; N=(2r+1)^2 except for boundary pixels.

mean_I_r = boxfilter(I(:, :, 1), r) ./ N;
mean_I_g = boxfilter(I(:, :, 2), r) ./ N;
mean_I_b = boxfilter(I(:, :, 3), r) ./ N;

% variance of I in each local patch: the matrix Sigma in Eqn (14).
% Note the variance in each local patch is a 3x3 symmetric matrix:
%           rr, rg, rb
%   Sigma = rg, gg, gb
%           rb, gb, bb
var_I_rr = boxfilter(I(:, :, 1).*I(:, :, 1), r) ./ N - mean_I_r .*  mean_I_r + eps; 
var_I_rg = boxfilter(I(:, :, 1).*I(:, :, 2), r) ./ N - mean_I_r .*  mean_I_g; 
var_I_rb = boxfilter(I(:, :, 1).*I(:, :, 3), r) ./ N - mean_I_r .*  mean_I_b; 
var_I_gg = boxfilter(I(:, :, 2).*I(:, :, 2), r) ./ N - mean_I_g .*  mean_I_g + eps; 
var_I_gb = boxfilter(I(:, :, 2).*I(:, :, 3), r) ./ N - mean_I_g .*  mean_I_b; 
var_I_bb = boxfilter(I(:, :, 3).*I(:, :, 3), r) ./ N - mean_I_b .*  mean_I_b + eps; 

inv_rr = var_I_gg .* var_I_bb - var_I_gb .^ 2;
inv_rg = var_I_gb .* var_I_rb - var_I_rg .* var_I_bb;
inv_rb = var_I_rg .* var_I_gb - var_I_gg .* var_I_rb;
inv_gg = var_I_rr .* var_I_bb - var_I_rb .^ 2;
inv_gb = var_I_rb .* var_I_rg - var_I_rr .* var_I_gb;
inv_bb = var_I_rr .* var_I_gg - var_I_rg .^ 2;
det_I = var_I_rr .* inv_rr + var_I_rg .* inv_rg + var_I_rb .* inv_rb;

inv_rr = inv_rr ./ det_I;
inv_rg = inv_rg ./ det_I;
inv_rb = inv_rb ./ det_I;
inv_gg = inv_gg ./ det_I;
inv_gb = inv_gb ./ det_I;
inv_bb = inv_bb ./ det_I;

q = zeros(hei, wid, nP);
ax = zeros(hei, wid, 3, nP);
bx = zeros(hei, wid, nP);
for iP = 1:nP,
    NW = boxfilter(ones(hei, wid) .* weight, r);

    mean_Ip_r = boxfilter(I(:, :, 1).*p(:,:,iP).*weight, r) ./ NW;
    mean_Ip_g = boxfilter(I(:, :, 2).*p(:,:,iP).*weight, r) ./ NW;
    mean_Ip_b = boxfilter(I(:, :, 3).*p(:,:,iP).*weight, r) ./ NW;
    mean_p = boxfilter(p(:,:,iP).*weight, r) ./ NW;

    % covariance of (I, p) in each local patch.
    cov_Ip_r = mean_Ip_r - mean_I_r .* mean_p;
    cov_Ip_g = mean_Ip_g - mean_I_g .* mean_p;
    cov_Ip_b = mean_Ip_b - mean_I_b .* mean_p;

    a_r = cov_Ip_r .* inv_rr + cov_Ip_g .* inv_rg + cov_Ip_b .* inv_rb;
    a_g = cov_Ip_r .* inv_rg + cov_Ip_g .* inv_gg + cov_Ip_b .* inv_gb;
    a_b = cov_Ip_r .* inv_rb + cov_Ip_g .* inv_gb + cov_Ip_b .* inv_bb;
    ax(:,:,1,iP) = a_r;
    ax(:,:,2,iP) = a_g;
    ax(:,:,3,iP) = a_b;

    b = mean_p - a_r .* mean_I_r - a_g .* mean_I_g - a_b .* mean_I_b; % Eqn. (15) in the paper;
    bx(:,:,iP) = b;

    q(:,:,iP) = (boxfilter(a_r, r).* I(:, :, 1)...
    + boxfilter(a_g, r).* I(:, :, 2)...
    + boxfilter(a_b, r).* I(:, :, 3)...
    + boxfilter(b, r)) ./ N;  % Eqn. (16) in the paper;
end
end