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

if nargin < 3
    arr_b = zeros(size(arr_y));
else
    arr_b = reshape(arr_b,1,[]);
end

% make h = 0 whenever y = 0 and vice versa
arr_h = bsxfun(@times,~isnan(arr_y),arr_h);
arr_y = bsxfun(@times,~isnan(arr_h),arr_y);

yh = nansum(bsxfun(@times,arr_y-arr_b,arr_h));
h2 = nansum(bsxfun(@power,arr_h,2));

c = yh / h2;

end % function

