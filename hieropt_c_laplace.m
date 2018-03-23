function [ c ] = hieropt_c_laplace( arr_y, arr_h, arr_noise, arr_b, c_mode )
% Computes the optimal proportionality parameter in a laplace model.
%
% Input:
%   arr_y     : observations
%   arr_h     : unscaled simulations
%   arr_noise : sigma2 noise
%   arr_b     : previously computed offset values
%   c_mode    : mode how c shall be computed. If 'absolute', just 1 is
%               returned, otherwise the optimal value computed.
%
% Output:
%   c         : optimal proportionality factor
%
% History:
%   2018/03/23: Yannik Schaelte

if strcmp(c_mode, 'absolute')
    c = 1;
    return;
end

% else compute optimal c

arr_y = reshape(arr_y,1,[]);
arr_h = reshape(arr_h,1,[]);
arr_b = reshape(arr_b,1,[]);
arr_recnoise = reshape(1./arr_noise,1,[]); % reciprocal noise

% unset values if not defined
bad_indices = ~isfinite(arr_y) | ~isfinite(arr_h) | ~isfinite(arr_recnoise) | ~isfinite(arr_b);
arr_y(bad_indices) = 0;
arr_h(bad_indices) = 0;
arr_b(bad_indices) = 0;
arr_recnoise(bad_indices) = 0; % here 1/inf = 0

end

