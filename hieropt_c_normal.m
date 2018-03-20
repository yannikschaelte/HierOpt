function [ c ] = hieropt_c_normal( arr_y, arr_h, arr_noise, arr_b, c_mode )
% [t,y,r,e]
% assumptions: resolution(b) = resolution(c), resolution(b) \supset
% resolution(sigma2)

if strcmp(c_mode,'absolute')
    c = 1;
    return;
end

% else compute optimal c

arr_y = reshape(arr_y,1,[]);
arr_h = reshape(arr_h,1,[]);
arr_recnoise = reshape(1./arr_noise,1,[]); % reciprocal noise

% unset values if not defined
bad_indices = ~isfinite(arr_y) | ~isfinite(arr_h) | ~isfinite(arr_recnoise);
arr_y(bad_indices) = 0;
arr_h(bad_indices) = 0;
arr_recnoise(bad_indices) = 0; % here 1/inf = 0

sum_yh = nansum((arr_y -arr_b) .* arr_h .* arr_recnoise);
sum_h2 = nansum(arr_h.^2 .*arr_recnoise);

c = sum_yh / sum_h2;

end % function

