function [ b ] = hieropt_b_normal( arr_y, arr_h, arr_noise, b_mode, c_mode )
% Computes the optimal offset parameter in a gaussian model.
%
% Input:
%   arr_y     : observations
%   arr_h     : unscaled simulations
%   arr_noise : sigma2 noise
%               we assume that either arr_noise = const., or 
%               noise_mode = 'absolute' (for the formulas to make sense)
%   b_mode    : mode how b shall be computed. If 'absolute', just 0 is
%               returned, otherwise the optimal value computed.
%   c_mode    : mode how c shall be computed, determines the formula for b
%
% Output:
%   b         : optimal offset value
%
% History:
%   2018/03/23: Yannik Schaelte

% [t,y,r,e]
% assumptions: resolution(b) = resolution(c), resolution(b) \supset
% resolution(sigma2), or sigma absolute (i.e. user-defined).

if strcmp(b_mode,'absolute')
    b = 0;
    return;
end

% else compute optimal b

arr_y = reshape(arr_y,1,[]);
arr_h = reshape(arr_h,1,[]);
arr_recnoise = reshape(1./arr_noise,1,[]); % reciprocal noise

% unset values if not defined
bad_indices = ~isfinite(arr_y) | ~isfinite(arr_h) | ~isfinite(arr_recnoise);
arr_y(bad_indices) = 0;
arr_h(bad_indices) = 0;
arr_recnoise(bad_indices) = 0; % here 1/inf = 0

sum_yh = nansum(arr_y .* arr_h .* arr_recnoise);
sum_h2 = nansum(arr_h.^2 .* arr_recnoise);
sum_h  = nansum(arr_h .* arr_recnoise);
sum_y  = nansum(arr_y .* arr_recnoise);
sum_recnoise = nansum(arr_recnoise);

% formula for b for known c

if strcmp(c_mode,'absolute')
    b = nansum((arr_y-arr_h).*arr_recnoise)/sum_recnoise;
    return;
end

% else formula for both b and c optimal

numerator = (arr_y - sum_yh*sum_h/sum_h2) / sum_recnoise;
denominator = 1 - (sum_h^2 / sum_h2) / sum_recnoise;

if abs(denominator) < eps
warning('hieropt:opt_b_normal: data not diverse enough to compute scalings this way. setting b = 0.');
    b = 0;
else
    b = numerator / denominator;
end

end

