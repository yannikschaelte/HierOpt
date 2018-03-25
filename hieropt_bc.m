function [ b, c] = hieropt_bc( arr_y, arr_h, arr_noise, b_mode, c_mode, distribution)
% Computes the optimal offset and proportionality parameters. If b_mode or
% c_mode equals 'absolute', the returned value is simply 0 or 1,
% respectively.
%
% Input:
%   arr_y        : observations
%   arr_h        : unscaled simulations
%   arr_noise    : sigma2 noise
%   b_mode       : mode how b shall be computed. If 'absolute', just 0 is
%                  returned, otherwise the optimal value computed.
%   c_mode       : mode how c shall be computed, determines the formula for
%                  b
%   distribution : 'normal','laplace'
%                  error model
% Output:
%   b         : optimal offset value
%   c         : optimal proportionality value
%
% History:
%   2018/03/24: Yannik Schaelte

switch distribution
    case 'normal'
        [b,c] = hieropt_bc_normal(arr_y,arr_h,arr_noise,b_mode,c_mode);
        
    case 'laplace'    
        [b,c] = hieropt_bc_laplace(arr_y,arr_h,arr_noise,b_mode,c_mode);
    
    otherwise
        error('hieropt:bc',"Distributon not recognized.");

end

