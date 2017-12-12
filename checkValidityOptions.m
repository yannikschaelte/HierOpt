function ret = checkValidityOptions(D,options)
% Function that checks whether the groups in the options are set properly. 
%
% Parameters:
%      D: see nlLH_fgh.m
%      options: see nlLH_fgh.m
%
% Return values:
%      ret: true if options are valid, false otherwise
%

n_e = size(D,2); %number of experiments
n_y = size(D(1).my,2); %number of observables
n_r = size(D(1).my,3); %number of replicates

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

%% CHECK VALIDITY OF OPTIONS
counter = 1;
mat_prop = nan(n_y,n_e);
for iy = 1:numel(options.obs_groups.proportionality)
    mat_prop(options.obs_groups.proportionality{iy},1) = counter;
    counter = counter+1;
end

counter = 0;
ystructure = mat_prop(:,1);
for ie = 1:numel(options.exp_groups.proportionality)
    for iie = 1:numel(options.exp_groups.proportionality{ie})
        mat_prop(:,options.exp_groups.proportionality{ie}(iie)) = ystructure + counter;
    end
    counter = counter + max(ystructure)+1;
end

counter = 1;
mat_var = nan(n_y,n_e);
for iy = 1:numel(options.obs_groups.variance)
    mat_var(options.obs_groups.variance{iy},1) = counter;
    counter = counter+1;
end

counter = 0;
ystructure = mat_var(:,1);
for ie = 1:numel(options.exp_groups.variance)
    for iie = 1:numel(options.exp_groups.variance{ie})
        mat_var(:,options.exp_groups.variance{ie}(iie)) = ystructure + counter;
    end
    counter = counter + max(ystructure)+1;
end
ret = true;
vals = unique(mat_prop);
for i = 1:numel(vals)
    if numel(unique(mat_var(mat_prop==vals(i)))) > 1
        ret = false;
    end
end