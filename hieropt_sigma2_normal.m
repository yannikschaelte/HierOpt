function [ sigma2 ] = hieropt_sigma2_normal( arr_y, arr_h, arr_c, arr_b )

arr_y = reshape(arr_y,1,[]);
arr_h = reshape(arr_h,1,[]);

if nargin < 3
    arr_c = ones(size(arr_y));
else
    arr_c = reshape(arr_c,1,[]);
end

if nargin < 4
    arr_b = ones(size(arr_y));
else
    arr_b = reshape(arr_b,1,[]);
end

% count how many observables (i.e. ~= nan) we have
count = sum(~isnan(arr_y));

% make h = 0 whenever y = 0 and vice versa
arr_h = bsxfun(@times,~isnan(arr_y),arr_h);
arr_y = bsxfun(@times,~isnan(arr_h),arr_y);

y_ch = bsxfun(@minus,arr_y,bsxfun(@times,arr_c,arr_h)+arr_b);

numerator = nansum(bsxfun(@power,y_ch,2));

sigma2 = numerator / count;

end

