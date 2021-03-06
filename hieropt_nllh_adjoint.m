function [ varargout ] = hieropt_nllh_adjoint( simfun, theta, D, amiOptions, scOptions )
% hieropt_nllh_adjoint uses the hierarchical approach to compute the nllh,
% and derivatives using adjoint sensitivity analysis.
% The objective function interface must be that of the DE integration 
% toolbox AMICI.
% The llh and derivatives returned by the simfun in field sim.llh must
% correspond to the error model specified in scOptions.distribution.
%
% Input:
%   simfun
%   theta
%   D
%   amiOptions
%   scOptions
%
% Output:
%   [nllh,snllh,s2nllh]
%
% History:
%   2018/01/12: Yannik Schaelte

if nargin <= 3
    amiOptions = amioption();
end
if nargin <= 4
    scOptions = struct();
end

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
    
    kappa_e = [D(ie).condition(:); zeros(n_obsGroups_notabs_b,1); ones(n_obsGroups_notabs_c,1)];
    
    sol = simfun(D(ie).t,theta,kappa_e,[],amiOptions);
    
    if (sol.status ~= 0)
        error('hieropt:nllh', "Could not integrate ODE.");
    end
    
    sim(ie).y = sol.y;
end

% optimal scalings
[b,c,noise,b_by_y,c_by_y,~] = hieropt_scalings(sim,D,scOptions);

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
    nllh = hieropt_nllh_nosim(sim,D,b,c,noise);
else
    % we need to perform the backward integration to obtain sensitivities
    
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
            
            sigma_re = sqrt(noise{ie}(:,:,ir));
%             sigma_re = reshape(sigma_re,1,[]);
%             sigma_re = repmat(sigma_re,length(D(ie).t),1);
            
            kappa_re = [D(ie).condition(:); b_re; c_re];
            
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
