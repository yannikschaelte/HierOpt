function [ b, c ] = hieropt_bc_laplace( arr_y, arr_h, arr_noise, b_mode, c_mode )
% Computes the optimal offset and proportionality parameters in a Laplace
% noise model. If b_mode or c_mode equals 'absolute', the returned value is 
% simply 0 or 1, respectively.
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

% else compute optimal value

arr_y = reshape(arr_y,1,[]);
arr_h = reshape(arr_h,1,[]);
arr_noise = reshape(arr_noise,1,[]);

dim = length(arr_y);

J = @(b,c) nansum( log(2*arr_noise) + abs( arr_y - (c*arr_h + b) ) ./ arr_noise);

if strcmp(b_mode, 'absolute')
   b = 0;
   
   if strcmp(c_mode, 'absolute')
       c = 1;
       return;
   end
   
   c_can = arr_y ./ arr_h;
   
   % determine minimum value
   J_can = nan(dim,1);
   for j = 1:dim
       if isfinite(c_can(j))
           J_can(j) = J(0,c_can(j));
       end
   end
   
   [~,i_min] = min(J_can);
   c = c_can(i_min);
end
    
end

