function [ noise ] = hieropt_noise( arr_y, arr_h, arr_b, arr_c, noise_mode, distribution)
% Computes the optimal noise parameter.
%
% Input:
%   arr_y        : observations
%   arr_h        : unscaled simulations
%   arr_b        : previously computed offset values
%   arr_c        : previously computed proportionality values
%   noise_mode   : 'absolute','single','multiple'
%                  noise mode
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
        noise = hieropt_noise_normal(arr_y,arr_h,arr_b,arr_c,noise_mode);
    case 'laplace'
        noise = hieropt_noise_laplace(arr_y,arr_h,arr_b,arr_c,noise_mode);
    otherwise
        error('hieropt:noise',"Distribution not recognized.");
        
end

