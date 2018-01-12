function [ b,c,sigma2,b_by_y,c_by_y,sigma2_by_y ] = hieropt_scalings(sim,data,scOptions)
% hieropt_scalings computes the optimal scalings (b,c,sigma2).
%
% Parameters:
%   scOptions:
%     .exp_groups
%       .bc_idxs
%       .sigma2_idxs
%     .obs_groups
%       .bc_idxs
%       .bc_mode       : 'multiple','single','absolute'
%       .sigma2_idxs   
%       .sigma2_mode   : 'multiple','single','absolute'
%
% Return Values:
%   b
%   c
%   sigma2
%   b_by_y
%   c_by_y
%   sigma2_by_y
% 
% History:
%   2018/01/12: Yannik Schaelte

% TODO we can also have user-defined sigma2s, giving different formulas for
% the other scalings

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

% compute numbers of not-absolute groups
n_obsGroups_notabs_b = 0;
n_obsGroups_notabs_c = 0;
n_obsGroups_notabs_sigma2 = 0;
for iyg = 1:n_obsGroups_bc
    if ~strcmp(scOptions.obs_groups.b_mode{iyg},'absolute')
        n_obsGroups_notabs_b = n_obsGroups_notabs_b + 1;
    end
    if ~strcmp(scOptions.obs_groups.c_mode{iyg},'absolute')
        n_obsGroups_notabs_c = n_obsGroups_notabs_c + 1;
    end
end
for iyg = 1:n_obsGroups_sigma2
    if ~strcmp(scOptions.obs_groups.sigma2_mode{iyg},'absolute')
        n_obsGroups_notabs_sigma2 + n_obsGroups_notabs_sigma2 + 1;
    end
end

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

    b_by_y{ie} = zeros(nt,n_obsGroups_notabs_b,nr);
    c_by_y{ie} = ones(nt,n_obsGroups_notabs_c,nr);
    sigma2_by_y{ie} = ones(nt,n_obsGroups_notabs_sigma2,nr);
end

%% OPTIMAL VALUES FOR THE PROPORTIONALITIES AND OFFSETS

for ieg = n_expGroups_bc
    
    ind_e = scOptions.exp_groups.bc_idxs{ieg};
    nr = size(data(ind_e(1)).my,3);
    
    i_obsGroups_notabs_b = 1;
    i_obsGroups_notabs_c = 1;
    
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
            tmp_b = hieropt_b_normal(arr_y,arr_h,b_mode,c_mode);
            for ie = ind_e
                b{ie}(:,ind_y,ind_r) = tmp_b;
                if ~strcmp(b_mode,'absolute')
                    b_by_y{ie}(:,i_obsGroups_notabs_b,ind_r) = tmp_b;
                    i_obsGroups_notabs_b = i_obsGroups_notabs_b + 1;
                end
            end
            arr_b = tmp_b*ones(size(arr_y));
            
            % compute optimal c
            tmp_c = hieropt_c_normal(arr_y,arr_h,arr_b,c_mode);
            for ie = ind_e
                c{ie}(:,ind_y,ind_r) = tmp_c;
                if ~strcmp(c_mode,'absolute')
                    c_by_y{ie}(:,i_obsGroups_notabs_c,ind_r) = tmp_c;
                    i_obsGroups_notabs_c = i_obsGroups_notabs_c + 1;
                end
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
            tmp_sigma2 = hieropt_sigma2_normal(arr_y,arr_h,arr_c,arr_b);
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
