function [ noise ] = hieropt_noise_normal( arr_y, arr_h, arr_b, arr_c )
% Computes noise=sigma2 in a Gaussian noise model.

arr_y = reshape(arr_y,1,[]);
arr_h = reshape(arr_h,1,[]);

% make h = 0 whenever y = 0 and vice versa
bad_indices = ~isfinite(arr_y) || ~isfinite(arr_h) || ~isfinite(arr_b) || ~isfinite(arr_c);
arr_h(bad_indides) = 0;
arr_y(bad_indices) = 0;
arr_b(bad_indices) = 0;
arr_c(bad_indices) = 0;
% count good indices
count = sum(~bad_indices);

y_ch = arr_y - (arr_c .* arr_h +arr_b);

numerator = nansum(y_ch.^2);

noise = numerator / count;

end

