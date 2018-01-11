function [ b,c,sigma2,b_by_y,c_by_y,sigma2_by_y ] = opt_scalings(sim,data,scOptions)
% INPUT
% scOptions:
%   .exp_groups
%     .bc_idxs
%     .sigma2_idxs
%   .obs_groups
%     .bc_idxs
%     .bc_mode       : 'multiple','single','absolute'
%     .sigma2_idxs   
%     .sigma2_mode   : 'multiple','single','absolute','user-specified'
%     .sigma2_values : used if sigma2_mode='user-specified', must be of
%                      fitting dimensions

%% PRELIMINARIES

% number of experiments
ne = size(data,2);
% number of observable (should be the same for all experiments, doesn't
% make sense really otherwise)
ny = size(data(1).my,2);

% fill unset fields with default values and perform sanity check
scOptions = sanityCheck(scOptions,ne,ny);

% initialization of b,c,sigma2
% return a single (c,sigma2) for every y,r,e (might return only one for
% every obs_group), since the adjoint computations will be possible only
% separatedly for each r,e

n_expGroups_bc = numel(scOptions.exp_groups.bc_idxs);
n_expGroups_sigma2 = numel(scOptions.exp_groups.sigma2_idxs);
n_obsGroups_bc     = numel(scOptions.obs_groups.bc_idxs);
n_obsGroups_sigma2 = numel(scOptions.obs_groups.sigma2_idxs);

b = cell(ne,1);
c = cell(ne,1);
sigma2 = cell(ne,1);

b_by_y = cell(ne,1);
c_by_y = cell(ne,1);
sigma2_by_y = cell(ne,1);

for ie = 1:ne
    nt = size(data(ie).my,1);
    nr = size(data(ie).my,3);
    
    b{ie} = zeros(nt,ny,nr);
    c{ie} = ones(nt,ny,nr);
    sigma2{ie} = ones(nt,ny,nr);

    b_by_y{ie} = zeros(n_obsGroups_bc,nr);
    c_by_y{ie} = ones(n_obsGroups_bc,nr);
    sigma2_by_y{ie} = ones(n_obsGroups_sigma2,nr);
end

%% OPTIMAL VALUES FOR THE PROPORTIONALITIES AND OFFSETS

for ieg = n_expGroups_bc
    
    ind_e = scOptions.exp_groups.bc_idxs{ieg};
    nr = size(data(ind_e(1)).my,3);
    
    for iyg = 1:n_obsGroups_bc
        
        ind_y = scOptions.obs_groups.bc_idxs{iyg};
        
        % create replicate groups
        % TODO we don't really need replicate groups, we can also have a
        % more generic approach at little cost
        b_mode = scOptions.obs_groups.b_mode{iyg};
        c_mode = scOptions.obs_groups.c_mode{iyg};
        if strcmp(b_mode,'multiple') && strcmp(c_mode,'single')...
            || strcmp(b_mode,'single') && strcmp(c_mode,'multiple')
            error('hieropt: not allowed combination of b_mode|c_mode single|multiple');
        end
        if strcmp(b_mode,'multiple') || strcmp(c_mode,'multiple')
            for ir=1:nr
                rep_groups.bc_idxs{ir}=ir;
            end
        elseif strcmp(b_mode,'single') || strcmp(c_mode,'single')
            rep_groups.bc_idxs{1}=1:nr;
        else
            % b,c have been set to 0,1 already
            continue;
        end
        n_repGroups_bc = numel(rep_groups.bc_idxs);
            
        for irg = 1:n_repGroups_bc
            
            ind_r = rep_groups.bc_idxs{irg};
            
            % create vectors with all entries
            arr_y = [];
            arr_h = [];
            for je = ind_e
               arr_y = [arr_y reshape(data(je).my(:,ind_y,ind_r),1,[])];
               arr_h = [arr_h reshape(sim(je).y(:,ind_y),1,[])];
            end
            
            % compute optimal b
            tmp_b = opt_b_normal(arr_y,arr_h,b_mode,c_mode);
            for ie = ind_e
                b{ie}(:,ind_y,ind_r) = tmp_b;
                b_by_y{ie}(:,iyg,ind_r) = tmp_b;
            end
            arr_b = tmp_b*ones(size(arr_y));
            
            % compute optimal c
            tmp_c = opt_c_normal(arr_y,arr_h,arr_b,c_mode);
            for ie = ind_e
                c{ie}(:,ind_y,ind_r) = tmp_c;
                c_by_y{ie}(:,iyg,ind_r) = tmp_c;
            end
        end
        
    end
end

%% OPTIMAL VALUES FOR THE VARIANCES
for ieg = 1:n_expGroups_sigma2
    
    ind_e = scOptions.exp_groups.sigma2_idxs{ieg};
    nr = size(data(ind_e(1)).my,3);
    
    for iyg = 1:n_obsGroups_sigma2
        
        ind_y = scOptions.obs_groups.sigma2_idxs{iyg};
        
        % create replicate groups
        % TODO we don't really need replicate groups, we can also have a
        % more generic approach at little cost
        switch scOptions.obs_groups.sigma2_mode{iyg}
            case 'multiple'
                for ir=1:nr
                    rep_groups.sigma2_idxs{ir}=ir;
                end
            case 'single'    
                rep_groups.sigma2_idxs{1}=1:nr;
            case 'absolute'
                % sigma2 has already been set to 1
                continue;
            case 'user-specified'
                error('TODO');
            otherwise
                error('could not resolve input');
        end
        n_repGroups_sigma2 = numel(rep_groups.sigma2_idxs);
        
        for irg = 1:n_repGroups_sigma2
            ind_r = rep_groups.sigma2_idxs{irg};
            
            % create vectors with all entries
            arr_y = [];
            arr_h = [];
            arr_c = [];
            arr_b = [];
            for ie = ind_e
               arr_y = [arr_y reshape(data(ie).my(:,ind_y,ind_r),1,[])];
               arr_h = [arr_h reshape(sim(ie).y(:,ind_y),1,[])];
               arr_c = [arr_c reshape(c{ie}(:,ind_y,ind_r),1,[])];
               arr_b = [arr_b reshape(b{ie}(:,ind_y,ind_r),1,[])];
            end
            
            % compute optimal sigma2
            tmp_sigma2 = opt_sigma2_normal(arr_y,arr_h,arr_c,arr_b);
            for ie = ind_e
                sigma2{ie}(:,ind_y,ind_r) = tmp_sigma2; 
                sigma2_by_y{ie}(:,iyg,ind_r) = tmp_sigma2;
            end
        end    
    end
end

end

function [ scOptions ] = sanityCheck(scOptions,ne,ny)

% default values for groupings: all together, absolute scalings,
% single sigma2
if ~isfield(scOptions,'exp_groups') || ~isfield(scOptions.exp_groups,'bc_idxs')
    scOptions.exp_groups.bc_idxs{1} = 1:ne;
end
if ~isfield(scOptions,'exp_groups') || ~isfield(scOptions.exp_groups,'sigma2_idxs')
    scOptions.exp_groups.sigma2_idxs{1} = 1:ne;
end
if ~isfield(scOptions,'obs_groups') || ~isfield(scOptions.obs_groups,'bc_idxs')
    scOptions.obs_groups.bc_idxs{1} = 1:ny;
    scOptions.obs_groups.b_mode{1} = 'absolute';
    scOptions.obs_groups.c_mode{1} = 'absolute';
end
if ~isfield(scOptions,'obs_groups') || ~isfield(scOptions.obs_groups,'sigma2_idxs')
    scOptions.obs_groups.sigma2_idxs{1} = 1:ny;
    scOptions.obs_groups.sigma2_mode{1} = 'absolute';
end

% check if modes of b,c are admissible

n_expGroups_bc = numel(scOptions.exp_groups.bc_idxs);
n_expGroups_sigma2 = numel(scOptions.exp_groups.sigma2_idxs);
n_obsGroups_bc     = numel(scOptions.obs_groups.bc_idxs);
n_obsGroups_sigma2 = numel(scOptions.obs_groups.sigma2_idxs);

for ieg = 1:n_expGroups_bc
    for iyg = 1:n_obsGroups_bc
        b_mode = scOptions.obs_groups.b_mode{iyg};
        c_mode = scOptions.obs_groups.c_mode{iyg};
        if strcmp(b_mode,'multiple') && strcmp(c_mode,'single')...
            || strcmp(b_mode,'single') && strcmp(c_mode,'multiple')
            error('hieropt: not allowed combination of b_mode|c_mode single|multiple');
        end
    end
end

% check if all values sharing b,c also share sigma2 (or sigma2
% not-optimized)

% TODO

% foreach bc_group check if exists! sigma2_group containing all indices

end


