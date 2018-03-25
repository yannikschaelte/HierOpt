function [ noise ] = hieropt_noise_laplace( arr_y, arr_h, arr_b, arr_c, noise_mode )
% Computes the optimal noise=sigma parameter in a laplacian model.
%
% Input:
%   arr_y      : observations
%   arr_h      : unscaled simulations
%   arr_b      : previously computed offset values
%   arr_c      : previously computed proportionality values
%   noise_mode : noise mode
%
% Output:
%   noise     : optimal noise=sigma2 value
%
% History:
%   2018/03/23: Yannik Schaelte

if strcmp(noise_mode,'absolute')
    noise = 1;
    return;
end

arr_y = reshape(arr_y,1,[]);
arr_h = reshape(arr_h,1,[]);
arr_b = reshape(arr_b,1,[]);
arr_c = reshape(arr_c,1,[]);

% make h = 0 whenever y = 0 and vice versa
bad_indices = ~isfinite(arr_y) | ~isfinite(arr_h) | ~isfinite(arr_b) | ~isfinite(arr_c);
arr_h(bad_indices) = 0;
arr_y(bad_indices) = 0;
arr_b(bad_indices) = 0;
arr_c(bad_indices) = 0;
% count good indices
count = sum(~bad_indices);

num = nansum( abs( arr_y - ( arr_c .* arr_h + arr_b ) ) );

noise = num / count;

end

