function [ b ] = hieropt_b_laplace( arr_y, arr_h, arr_noise, b_mode, c_mode )
% Computes the optimal offset parameter in a gaussian model.
%
% Input:
%   arr_y     : observations
%   arr_h     : unscaled simulations
%   arr_noise : sigma2 noise
%   b_mode    : mode how b shall be computed. If 'absolute', just 0 is
%               returned, otherwise the optimal value computed.
%   c_mode    : mode how c shall be computed, determines the formula for b
%
% Output:
%   b         : optimal offset value
%
% History:
%   2018/03/23: Yannik Schaelte

b = 0;

end

