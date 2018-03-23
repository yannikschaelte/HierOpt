function [ noise ] = hieropt_noise( arr_y, arr_h, arr_b, arr_c, distribution)
% Computes the optimal noise parameter.
%
% Input:
%   arr_y        : observations
%   arr_h        : unscaled simulations
%   arr_b        : previously computed offset values
%   arr_c        : previously computed proportionality values
%   distribution : 'normal','laplace'
%                  error model
%
% Output:
%   noise     : optimal noise value
%
% History:
%   2018/03/23: Yannik Schaelte

switch distribution
    case 'normal'
        noise = hieropt_noise_normal(arr_y,arr_h,arr_b,arr_c);
    case 'laplace'
        noise = hieropt_noise_laplace(arr_y,arr_h,arr_b,arr_c);
    otherwise
        error('Distribution not recognized.');
        
end

