function [ b ] = hieropt_b_normal( arr_y, arr_h, arr_noise, b_mode, c_mode )
% [t,y,r,e]
% assumptions: resolution(b) = resolution(c), resolution(b) \supset
% resolution(sigma2)

if strcmp(b_mode,'absolute')
    b = 0;
    return;
end

% else compute optimal b

arr_y = reshape(arr_y,1,[]);
arr_h = reshape(arr_h,1,[]);

% count how many observables (i.e. ~= nan) we have
count = sum(~isnan(arr_y));

% make h = 0 whenever y = 0 and vice versa
arr_h = bsxfun(@times,~isnan(arr_y),arr_h);
arr_y = bsxfun(@times,~isnan(arr_h),arr_y);

yh = nansum(bsxfun(@times,arr_y,arr_h));
h2 = nansum(bsxfun(@power,arr_h,2));

if strcmp(c_mode,'absolute')
    b = nansum(arr_y-arr_h)/count;
    return;
end

% else both b and c optimal

numerator = nansum(arr_y-(yh/h2)*arr_h)/count;
denominator = 1 - nansum((nansum(arr_h)/h2)*arr_h)/count;

if abs(denominator) < eps
warning('hieropt:opt_b_normal: data not diverse enough to compute scalings this way');
    b = 0;
else
    b = numerator / denominator;
end

end

