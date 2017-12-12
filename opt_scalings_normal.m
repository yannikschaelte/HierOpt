function [ c,sigma2,c_by_y,sigma2_by_y,b ] = opt_scalings_normal(sim,data,options)
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
n_t = size(data(1).my,1);
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

distr = 'normal';

% initialization of proportionality and variance parameters
% we return a single (c,sigma2) for every y,r,e (might return only one for
% every obs_group), since the adjoint computations will be possible only
% separatedly for each r,e
c = zeros(n_t,n_y,n_r,n_e); % vector including proportionality factors
b = zeros(n_t,n_y,n_r,n_e);
sigma2 = zeros(n_t,n_y,n_r,n_e); % vector including variances

% group scalings by observables
nObsGroupsProportionality = numel(options.obs_groups.proportionality);
nObsGroupsVariance        = numel(options.obs_groups.variance);
c_by_y = zeros(1,nObsGroupsProportionality,n_r,n_e); % vector including proportionality factors
b_by_y = zeros(1,nObsGroupsProportionality,n_r,n_e);
sigma2_by_y = zeros(1,nObsGroupsVariance,n_r,n_e); % vector including variances

%% OPTIMAL VALUES FOR THE PROPORTIONALITY FACTORS AND OFFSETS
for ie = 1:numel(options.exp_groups.proportionality)
    
    ind_e = options.exp_groups.proportionality{ie};
    n_r = size(data(ind_e(1)).my,3);
    
    for iy = 1:numel(options.obs_groups.proportionality)
        ind_y = options.obs_groups.proportionality{iy};
        
        % create replicate groups
        % TODO we don't really need replicate groups, we can also have a
        % more generic approach at little cost
        switch options.obs(ind_y(1)).proportionality
            case 'multiple'
                for ir=1:n_r
                    rep_groups.proportionality{ir}=ir;
                end
            case 'single'    
                rep_groups.proportionality{1}=1:n_r;
            case 'absolute'
                c(:,ind_y,:,ind_e) = 1;
                c_by_y(:,iy,:,ind_e) = 1;
                continue;
            otherwise
                error('could not resolve input');
        end
            
        for ir = 1:numel(rep_groups.proportionality)
            ind_r = rep_groups.proportionality{ir};
            
            % create vectors with all entries
            arr_y = [];
            arr_h = [];
            for je = ind_e
               arr_y = [arr_y reshape(data(je).my(:,ind_y,ind_r),1,[])];
               arr_h = [arr_h reshape(sim(je).y(:,ind_y),1,[])];
            end
            
            b(:,ind_y,ind_r,ind_e) = opt_c_normal(arr_y,arr_h);
            arr_b = zeros(size(arr_y));
            
            % compute optimal cs
            tempc = opt_c_normal(arr_y,arr_h,arr_b);
            c(:,ind_y,ind_r,ind_e) = tempc; 
            c_by_y(:,iy,ind_r,ind_e) = tempc;
        end
        
    end
end

% c = c(1,:,:,:);

%% OPTIMAL VALUES FOR THE VARIANCES
for ie = 1:numel(options.exp_groups.variance)
    
    ind_e = options.exp_groups.variance{ie};
    n_r = size(data(ind_e(1)).my,3);
    
    for iy = 1:numel(options.obs_groups.variance)
        ind_y = options.obs_groups.variance{iy};
        
        % create replicate groups
        % TODO we don't really need replicate groups, we can also have a
        % more generic approach at little cost
        switch options.obs(ind_y(1)).variance
            case 'multiple'
                for ir=1:n_r
                    rep_groups.variance{ir}=ir;
                end
            case 'single'    
                rep_groups.variance{1}=1:n_r;
            case 'absolute'
                sigma2(:,ind_y,:,ind_e) = 1;
                sigma2_by_y(:,iy,:,ind_e) = 1;
                continue;
            otherwise
                error('could not resolve input');
        end
            
        for ir = 1:numel(rep_groups.variance)
            ind_r = rep_groups.variance{ir};
            
            % create vectors with all entries
            arr_y = [];
            arr_h = [];
            arr_c = [];
            for je = ind_e
               arr_y = [arr_y reshape(data(je).my(:,ind_y,ind_r),1,[])];
               arr_h = [arr_h reshape(sim(je).y(:,ind_y),1,[])];
               arr_c = [arr_c reshape(c(:,ind_y,ind_r,je),1,[])];
            end
            
%             arr_c = c(:,ind_y,ind_r,ind_e);
            arr_c = reshape(arr_c,1,[]);
            
            % compute optimal cs
            tempsigma2 = opt_sigma2_normal(arr_y,arr_h,arr_c,arr_b);
            sigma2(:,ind_y,ind_r,ind_e) = tempsigma2; 
            sigma2_by_y(:,iy,ind_r,ind_e) = tempsigma2;
        end
        
    end
    
    
%     for iy = 1:numel(options.obs_groups.variance)
%         ind_y = options.obs_groups.variance{iy};
%         ind_e = options.exp_groups.variance{ie};
%         tempsigma2 = ...
%             optim_sigma2(ind_y,sim(ind_e),data(ind_e),distr,options,c(:,:,:,ind_e));
%         for ir = 1:n_r
%             sigma2(:,ind_y,ir,ind_e) = tempsigma2(:,:,ir);
%             sigma2_by_y(:,iy,ir,ind_e) = tempsigma2(:,:,ir);
%         end
%         
%     end
end

% for the further processing, only return one c, since c will be the same
% for all time points, if we for example want to use it in the standard way
% with amici and adjoints
b = b(1,:,:,:);
c = c(1,:,:,:);
% sigma2 = sigma2(1,:,:,:);
