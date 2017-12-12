function [ c,sigma2,c_by_y,sigma2_by_y ] = scalings(sim,data,options)
% This function computes the value of the negative-log-likelihood function,
% its gradient and its Hessian at theta
%
% USAGE:
% [...]                     = nlLH_fgh(simulation,D,'normal')
% [...]                     = nlLH_fgh(simulation,D,'normal',options,1)
% [nllH]                    = nlLH_fgh(...)
% [nllH,gradnlLH]           = nlLH_fgh(...)
% [nllH,gradnlLH,FIMnlLH]   = nlLH_fgh(...)
%
% Parameters:
%   simulation: (1 x #experiments) struct with fields
%       * y: simulation of the output for the different experiments for a theta
%            in which the values nlLH, gradnlLH, FIMnlLH will be computed
%       * sy: simulation of sy, the sensitivities of the output vector, for the different experiments
%            for a theta in which the values nlLH, gradnlLH, FIMnlLH will be computed
%   D: (1 x #experiments) struct containing data with two fields
%       * t: time points
%       * my: number time points x number observables x number replicates
%   distr: indicates distribution of measurement noise,
%         = 'normal'
%         = 'laplace'
%         = 'log-normal'
%         = 'log-laplace'
%   options.exp_groups:
%               * variance: struct with each struct containing the indices
%                   of the experiments that share a variance parameter
%               * proportionality: struct with each struct containing the indices
%                   of the experiments that share a proportionality parameter
%   options.obs_groups:
%               * variance: struct with each struct containing the indices
%                   of the observables that share a variance parameter
%               * proportionality: struct with each struct containing the indices
%                   of the observables that share a proportionality parameter
%    (the defaults for the above options are that all parameters are shared
%    across observables and experiments)
%   options.obs: (1 x #observables) struct with fields
%               * variance:
%                           = 'single': the variance is computed for all replicate together
%                           = 'multiple': the variance is computed for each replicate separately
%               * proportionality:
%                           = 'single': proportionality factor is computed
%                                   for all replicates together
%                           = 'multiple': proportionality factor is
%                                   computed for each replicate separately
%                           = 'absolute': proportionality factor is 1
%   save_pv: indicates whether the proportionality factors and variances are
%            saved in a .mat file (1) or not (0)
%
% Return values:
%   nlLh: value of the negative-log-likelihood function in theta
%   gradnlLH: value of gradient of the negative-log-likelihood function
%                 in theta
%   FIMnlLH: approximation of the Hessian of the negative-log-likelihood function in theta
%            using the FIM
%
% Note: all oberservables/experiments that share a proportionality parameter
%       need to also share the variance parameter!
%
%% INITIALIZATION

n_e = size(data,2); %number of experiments
n_y = size(data(1).my,2); %number of observables
n_r = size(data(1).my,3); %number of replicates

% default values for groupings: all together
if ~isfield(options,'exp_groups')
    options.exp_groups.proportionality{1} = [1:n_e];
    options.exp_groups.variance{1} = [1:n_e];
else
    if ~isfield(options.exp_groups,'proportionality')
        options.exp_groups.proportionality{1} = [1:n_e];
    end
    if ~isfield(options.exp_groups,'variance')
        options.exp_groups.variance{1} = [1:n_e];
    end
end
if ~isfield(options,'obs_groups')
    options.obs_groups.proportionality{1} = [1:n_y];
    options.obs_groups.variance{1} = [1:n_y];
else
    if ~isfield(options.obs_groups,'proportionality')
        options.obs_groups.proportionality{1} = [1:n_y];
    end
    if ~isfield(options.obs_groups,'variance')
        options.obs_groups.variance{1} = [1:n_y];
    end
end

nlLH = 0;

distr = 'normal';

% initialization of proportionality and variance parameters
% we return a single (c,sigma2) for every y,r,e (might return only one for
% every obs_group), since the adjoint computations will be possible only
% separatedly for each r,e
c = zeros(1,n_y,n_r,n_e); % vector including proportionality factors
sigma2 = zeros(1,n_y,n_r,n_e); % vector including variances

% group scalings by observables
nObsGroupsProportionality = numel(options.obs_groups.proportionality);
nObsGroupsVariance        = numel(options.obs_groups.variance);
c_by_y = zeros(1,nObsGroupsProportionality,n_r,n_e); % vector including proportionality factors
sigma2_by_y = zeros(1,nObsGroupsVariance,n_r,n_e); % vector including variances

%% OPTIMAL VALUES FOR THE PROPORTIONALITY FACTORS
for ie = 1:numel(options.exp_groups.proportionality)
    for iy = 1:numel(options.obs_groups.proportionality)
        ind_y = options.obs_groups.proportionality{iy};
        ind_e = options.exp_groups.proportionality{ie};
        
        tempc = optim_c(ind_y,sim(ind_e),data(ind_e),distr,options);
        for ir = 1:n_r
            c(:,ind_y,ir,ind_e) = tempc(:,:,ir);
            c_by_y(:,iy,ir,ind_e) = tempc(:,:,ir);
        end        
    end
end

%% OPTIMAL VALUES FOR THE VARIANCES
for ie = 1:numel(options.exp_groups.variance)
    for iy = 1:numel(options.obs_groups.variance)
        ind_y = options.obs_groups.variance{iy};
        ind_e = options.exp_groups.variance{ie};
        tempsigma2 = ...
            optim_sigma2(ind_y,sim(ind_e),data(ind_e),distr,options,c(:,:,:,ind_e));
        for ir = 1:n_r
            sigma2(:,ind_y,ir,ind_e) = tempsigma2(:,:,ir);
            sigma2_by_y(:,iy,ir,ind_e) = tempsigma2(:,:,ir);
        end
        
    end
end