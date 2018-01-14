function [ varargout ] = hieropt_nllh_adjoint( simfun,theta,D,amiOptions,scOptions )
% hieropt_nllh_adjoint uses the hierarchical approach to compute the nllh.
% Derivatives are computed using the adjoint approach.
% For an example usage see
% examples/HierOpt_Examples/jakstat_small/nllh_jakstat_hierarchical_adjoint
% (_offsets).
%
% Parameters:
%   simfun
%   theta
%   D
%   amiOptions
%   scOptions
%
% Return Values:
%   [nllh,snllh,s2nllh]
%
% History:
%   2018/01/12: Yannik Schaelte

ne = size(D,2);
n_obsGroups_bc     = numel(scOptions.obs_groups.bc_idxs);

% compute numbers of not-absolute groups
n_obsGroups_notabs_b = 0;
n_obsGroups_notabs_c = 0;
for iyg = 1:n_obsGroups_bc
    if ~strcmp(scOptions.obs_groups.b_mode{iyg},'absolute')
        n_obsGroups_notabs_b = n_obsGroups_notabs_b + 1;
    end
    if ~strcmp(scOptions.obs_groups.c_mode{iyg},'absolute')
        n_obsGroups_notabs_c = n_obsGroups_notabs_c + 1;
    end
end

% forward simulation

amiOptions.sensi = 0;
sim = struct([]);
for ie = 1:ne
    
    kappa_e = [D(ie).k(:); zeros(n_obsGroups_notabs_b,1); ones(n_obsGroups_notabs_c,1)];
    
    sol = simfun(D(ie).t,theta,kappa_e,[],amiOptions);
    
    if (sol.status ~= 0)
        error('hieropt: could not integrate ODE.');
    end
    
    sim(ie).y = sol.y;
    if nargout > 1
        sim(ie).sy = sol.sy;
    end
end

% optimal scalings
[b,c,sigma2,b_by_y,c_by_y,sigma2_by_y] = hieropt_scalings(sim,D,scOptions);

% nllh, snllh, s2nllh

% set correct amiOptions
amiOptions.sensi_meth = 'adjoint';
switch nargout
    case 1
        amiOptions.sensi = 0;
    case 2
        amiOptions.sensi = 1;
    case 3
        amiOptions.sensi = 2;
end

if nargout == 1
    nllh = hieropt_nllh_forward(false,sim,D,b,c,sigma2);
else
    nllh = 0;
    n_theta = length(theta);
    snllh = zeros(n_theta,1);
    if nargout > 2
        s2nllh = zeros(n_theta);
    end
    
    for ie = 1:ne
        nr = size(D(ie).Y,3);
        for ir = 1:nr
            clear amiData
            
            b_re = b_by_y{ie}(1,:,ir);
            b_re = reshape(b_re,[],1);
            
            c_re = c_by_y{ie}(1,:,ir);
            c_re = reshape(c_re,[],1);
            
            sigma_re = sqrt(sigma2{ie}(1,:,ir));
            sigma_re = reshape(sigma_re,1,[]);
            sigma_re = repmat(sigma_re,length(D(ie).t),1);
            
            kappa_re = [D(ie).k(:); b_re; c_re];
            
            amiData.t = D(ie).t;
            amiData.Y = D(ie).Y(:,:,ir);
            amiData.Sigma_Y = sigma_re;
            amiData.condition = kappa_re;
            amiData = amidata(amiData);
            
            sol = simfun(D(ie).t,theta,kappa_re,amiData,amiOptions);
            
            if sol.status ~= 0
                error('hieropt: could not integrate ODE.');
            end
            
            nllh = nllh - sol.llh;
            snllh = snllh - sol.sllh;
            if nargout > 2
                s2nllh = s2nllh - sol.s2llh;
            end
            
        end
    end
end

varargout{1} = nllh;
if nargout > 1
    varargout{2} = snllh;
    if nargout > 2
        varargout{3} = s2nllh;
    end
end

end
