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

% else compute optimal value

arr_y = reshape(arr_y,1,[]);
arr_h = reshape(arr_h,1,[]);
arr_noise = reshape(arr_noise,1,[]);

dim = length(arr_y);

J = @(b,c) nansum( log(2*arr_noise) + abs( arr_y - (c*arr_h + b) ) ./ arr_noise);

if strcmp(b_mode, 'absolute')
    b = 0;
    
    if strcmp(c_mode, 'absolute')
        c = 1;
        return;
    end
    
    c_can = (arr_y-b) ./ arr_h;
    
    % determine minimum value
    J_can = nan(dim,1);
    for j = 1:dim
        if isfinite(c_can(j))
            J_can(j) = J(b,c_can(j));
        end
    end
    [~,i_min] = min(J_can);
    c = c_can(i_min);
    
    if ~isfinite(c)
        % means arr_h == 0 | arr_y, arr_h not finite anyway
        c = 1;
    end
    
else
    if strcmp(c_mode, 'absolute')
        c = 1;
        
        b_can = arr_y - c*arr_h;
        J_can = nan(dim,1);
        for j = 1:dim
            J_can(j) = J(b_can(j),c);
        end
        [~,i_min] = min(J_can);
        b = b_can(j);
        
    else
        % both b, c to be optimized
        J_min = inf;
        b_min = nan;
        c_min = nan;
        
        for ic = 1:dim
            for ib = 1:dim
                factor = 1 - arr_h(ib)/arr_h(ic);
                if ~isfinite(factor) || abs(factor) < eps
                    b_tmp = 0;
                else
                    b_tmp = (arr_y(ib) - arr_h(ib)/arr_h(ic)*arr_y(ic)) / factor;
                end
                c_tmp = (arr_y(ic) - b_tmp) / arr_h(ic);
                
                % check J
                J_tmp = J(b_tmp,c_tmp);
                if J_tmp < J_min
                    J_min = J_tmp;
                    b_min = b_tmp;
                    c_min = c_tmp;
                end
            end
        end
        
        b = b_min;
        c = c_min;
        
        if ~isfinite(b)
            b = 0;
        end
        
        if ~isfinite(c)
            c = 1;
        end
        
    end
    
end

