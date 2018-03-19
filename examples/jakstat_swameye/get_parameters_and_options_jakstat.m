function [parameters,options] = getParametersAndOptions_jakstat(approach)

n_starts = 500;

options.MS = PestoOptions();
options.MS.n_starts = n_starts; % actually 500
options.MS.mode = 'text';
options.MS.localOptimizer = 'fmincon';
options.MS.localOptimizerOptions = optimset('algorithm','interior-point',...
    'display','iter',...
    'GradObj','on',...
    'MaxIter',2000,...
    'TolFun',1e-10,...
    'TolX',1e-10,...
    'MaxFunEvals',2000,...
    'PrecondBandWidth', inf);
options.MS.obj_type = 'negative log-posterior';

n_par = 17; % maximum number of parameters
min_par = -5*ones(n_par,1);
max_par = 3*ones(n_par,1);
max_par(4)  =  6;
max_par(2)  =  6;
min_par(10) = -6;
min_par(4)  = -3;
min_par(2)  = -3;

par0 = bsxfun(@plus,min_par,bsxfun(@times,max_par - min_par, lhsdesign(n_starts,n_par,'smooth','off')'));
   
switch approach
    case 'standard'
        n_par = 17;
        
    case 'hierarchical-adjoint'
        n_par = 12;
        
        sc.exp_groups.bc_idxs = {1};
        sc.exp_groups.noise_idxs = {1};
        
        sc.obs_groups.bc_idxs = {1,2,3};
        sc.obs_groups.b_mode = {'absolute','absolute','absolute'};
        sc.obs_groups.c_mode = {'multiple','multiple','absolute'};
        sc.obs_groups.noise_idxs = {1,2,3};
        sc.obs_groups.noise_mode = {'multiple','multiple','multiple'};
        
        sc.distribution = 'normal';
        
        options.sc = sc;
        
    case 'hierarchical-forward'
        n_par = 12;
        
        sc.exp_groups.bc_idxs = {1};
        sc.exp_groups.noise_idxs = {1};
        
        sc.obs_groups.bc_idxs = {1,2,3};
        sc.obs_groups.b_mode = {'absolute','absolute','absolute'};
        sc.obs_groups.c_mode = {'multiple','multiple','absolute'};
        sc.obs_groups.noise_idxs = {1,2,3};
        sc.obs_groups.noise_mode = {'multiple','multiple','multiple'};
        
        sc.distribution = 'normal';
        
        options.sc = sc;
        
    case 'hierarchical-forward-offsets'
        n_par = 10;

        sc.exp_groups.bc_idxs = {1};
        sc.exp_groups.noise_idxs = {1};
        
        sc.obs_groups.bc_idxs = {1,2,3};
        sc.obs_groups.b_mode = {'multiple','multiple','absolute'};
        sc.obs_groups.c_mode = {'multiple','multiple','absolute'};
        sc.obs_groups.noise_idxs = {1,2,3};
        sc.obs_groups.noise_mode = {'multiple','multiple','multiple'};
        
        sc.distribution = 'normal';
        
        options.sc = sc;
        
    case 'hierarchical-adjoint-offsets'
        n_par = 10;
        
        
        sc.exp_groups.bc_idxs = {1};
        sc.exp_groups.noise_idxs = {1};
        
        sc.obs_groups.bc_idxs = {1,2,3};
        sc.obs_groups.b_mode = {'multiple','multiple','absolute'};
        sc.obs_groups.c_mode = {'multiple','multiple','absolute'};
        sc.obs_groups.noise_idxs = {1,2,3};
        sc.obs_groups.noise_mode = {'multiple','multiple','multiple'};
        
        sc.distribution = 'normal';
        
        options.sc = sc;
        
end

parameters.number = n_par;
parameters.min = min_par(1:n_par,1);
parameters.max = max_par(1:n_par,1);
parameters.guess = par0(1:n_par,1:n_starts);

end
