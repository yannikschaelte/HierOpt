function [ c ] = hieropt_c( arr_y, arr_h, arr_noise, arr_b, c_mode, distribution)
% Computes the optimal proportionality parameter.
%
% Input:
%   arr_y        : observations
%   arr_h        : unscaled simulations
%   arr_noise    : noise
%   arr_b        : previously computed offset values
%   c_mode       : mode how c shall be computed. If 'absolute', just 1 is
%                  returned, otherwise the optimal value computed.
%   distribution : 'normal','laplace'
%                  error model
%
% Output:
%   c         : optimal proportionality factor
%
% History:
%   2018/03/23: Yannik Schaelte

switch distribution
    case 'normal'
        c = hieropt_c_normal(arr_y,arr_h,arr_noise,arr_b,c_mode);
    case 'laplace'
        c = hieropt_c_laplace(arr_y,arr_h,arr_noise,arr_b,c_mode);
    otherwise
        error('Distribution invalid.');
        
end

