function [ b, c ] = hieropt_bc_normal( arr_y, arr_h, arr_noise, b_mode, c_mode )
% Computes the optimal offset and proportionality parameters in a Gaussian
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

b = hieropt_b_normal(arr_y,arr_h,arr_noise,b_mode,c_mode);

arr_b = b*ones(size(arr_y));

c = hieropt_c_normal(arr_y,arr_h,arr_noise,arr_b,c_mode);

end

