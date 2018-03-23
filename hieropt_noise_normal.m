function [ noise ] = hieropt_noise_normal( arr_y, arr_h, arr_b, arr_c )
% Computes the optimal noise=sigma2 parameter in a gaussian model.
%
% Input:
%   arr_y     : observations
%   arr_h     : unscaled simulations
%   arr_b     : previously computed offset values
%   arr_c     : previously computed proportionality values
%
% Output:
%   noise     : optimal noise=sigma2 value
%
% History:
%   2018/03/23: Yannik Schaelte

arr_y = reshape(arr_y,1,[]);
arr_h = reshape(arr_h,1,[]);

% make h = 0 whenever y = 0 and vice versa
bad_indices = ~isfinite(arr_y) | ~isfinite(arr_h) | ~isfinite(arr_b) | ~isfinite(arr_c);
arr_h(bad_indices) = 0;
arr_y(bad_indices) = 0;
arr_b(bad_indices) = 0;
arr_c(bad_indices) = 0;
% count good indices
count = sum(~bad_indices);

y_ch = arr_y - (arr_c .* arr_h +arr_b);

numerator = nansum(y_ch.^2);

noise = numerator / count;

end

