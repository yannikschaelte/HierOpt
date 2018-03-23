function [ b ] = hieropt_b_laplace( arr_y, arr_h, arr_noise, b_mode, c_mode )
% Computes the optimal offset parameter in a gaussian model.
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

if strcmp(b_mode, 'absolute')
    b = 0;
    return;
end

% else compute optimal value

arr_y = reshape(arr_y,1,[]);
arr_h = reshape(arr_h,1,[]);
arr_recnoise = reshape(1./arr_noise,1,[]); % reciprocal noise

% unset values if not defined
bad_indices = ~isfinite(arr_y) | ~isfinite(arr_h) | ~isfinite(arr_recnoise);
arr_y(bad_indices) = 0;
arr_h(bad_indices) = 0;
arr_recnoise(bad_indices) = 0; % here 1/inf = 0

dJdb = @(b) - sum

if strcmp(c_mode, 'absolute')
    b_can = arr_y - arr_h;
    
end
    
end

