function [ b ] = opt_b_normal( arr_y, arr_h )
% [t,y,r,e]
% assumptions: resolution(b) = resolution(c), resolution(b) \supset
% resolution(sigma2)

arr_y = reshape(arr_y,1,[]);
arr_h = reshape(arr_h,1,[]);

% count how many observables (i.e. ~= nan) we have
count = sum(~isnan(arr_y));

% make h = 0 whenever y = 0 and vice versa
arr_h = bsxfun(@times,~isnan(arr_y),arr_h);
arr_y = bsxfun(@times,~isnan(arr_h),arr_y);

yh = nansum(bsxfun(@times,arr_y,arr_h));
h2 = nansum(bsxfun(@power,arr_h,2));

numerator = nansum(y-(yh/h2)*h)/count;
denominator = 1 - nansum((nansum(h)/h2)*h)/count;

if abs(denominator) < eps
    warning('hieropt:opt_b_normal', 'data not diverse enough to compute scalings this way');
end

b = numerator / denominator;

end

