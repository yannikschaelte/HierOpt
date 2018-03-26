function [ c ] = hieropt_c_normal( arr_y, arr_h, arr_noise, arr_b, c_mode )
% Computes the optimal proportionality parameter in a gaussian model.
%
% Input:
%   arr_y     : observations
%   arr_h     : unscaled simulations
%   arr_noise : sigma2 noise
%               we assume that either arr_noise = const., or 
%               noise_mode = 'absolute' (for the formulas to make sense)
%   arr_b     : previously computed offset values
%   c_mode    : mode how c shall be computed. If 'absolute', just 1 is
%               returned, otherwise the optimal value computed.
%
% Output:
%   c         : optimal proportionality factor
%
% History:
%   2018/03/23: Yannik Schaelte

% [t,y,r,e]
% assumptions: resolution(b) = resolution(c), resolution(b) \supset
% resolution(sigma2)

if strcmp(c_mode, 'absolute')
    c = 1;
    return;
end

% else compute optimal c

arr_y = reshape(arr_y,1,[]);
arr_h = reshape(arr_h,1,[]);
arr_b = reshape(arr_b,1,[]);
arr_recnoise = reshape(1./arr_noise,1,[]); % reciprocal noise

% unset values if not defined
bad_indices = ~isfinite(arr_y) | ~isfinite(arr_h) | ~isfinite(arr_recnoise) | ~isfinite(arr_b);
arr_y(bad_indices) = 0;
arr_h(bad_indices) = 0;
arr_b(bad_indices) = 0;
arr_recnoise(bad_indices) = 0; % here 1/inf = 0

sum_yh = nansum((arr_y - arr_b) .* arr_h .* arr_recnoise);
sum_h2 = nansum(arr_h.^2 .*arr_recnoise);

% if abs(sum_h2) < eps
%     c = 1;
%     return;
% end

c = sum_yh / sum_h2;

end % function

