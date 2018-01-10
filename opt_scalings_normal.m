function [ c,sigma2,c_by_y,sigma2_by_y,b,b_by_y ] = opt_scalings_normal(sim,data,options)
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
                b(:,ind_y,:,ind_e) = 0;
                b_by_y(:,iy,:,ind_e) = 0;
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
            
            tempb = 0;%opt_b_normal(arr_y,arr_h);
            b(:,ind_y,ind_r,ind_e) = tempb;
            b_by_y(:,iy,ind_r,ind_e) = tempb;
%             arr_b = zeros(size(arr_y));
            arr_b = tempb*ones(size(arr_y));
            
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
            arr_b = [];
            for je = ind_e
               arr_y = [arr_y reshape(data(je).my(:,ind_y,ind_r),1,[])];
               arr_h = [arr_h reshape(sim(je).y(:,ind_y),1,[])];
               arr_c = [arr_c reshape(c(:,ind_y,ind_r,je),1,[])];
               arr_b = [arr_b reshape(b(:,ind_y,ind_r,je),1,[])];
            end
            
%             arr_c = c(:,ind_y,ind_r,ind_e);
%             arr_c = reshape(arr_c,1,[]);
% arr_b = zeros(size(arr_c));
            
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
