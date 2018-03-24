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

f_dJdb = @(b, c) - nansum(arr_recnoise * sign(arr_y - c*arr_h - b));
dim = length(arr_y);

if strcmp(c_mode, 'absolute')
    b_can = arr_y - arr_h;
    dJdB_can = zeros(dim+1,1);
    for j = 1:dim+1
        dJdB_can(j) = f_dJdb();%TODO
    end
end
    
end

