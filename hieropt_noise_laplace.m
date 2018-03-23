function [ noise ] = hieropt_noise_laplace( arr_y, arr_h, arr_b, arr_c )
% Computes the optimal noise=sigma parameter in a laplacian model.
%
% Input:
%   arr_y     : observations
%   arr_h     : unscaled simulations
%   arr_b     : previously computed offset values
%   arr_c     : previously computed proportionality values
%
% Output:
%   noise     : optimal noise=sigma2 value
%
% History:
%   2018/03/23: Yannik Schaelte

noise = 1;

end

