function [parameters,options] = get_parameters_and_options_erbb()

n_starts = 10;

options.MS = PestoOptions();
options.MS.n_starts = n_starts;
options.MS.mode = 'text';
options.MS.localOptimizer = 'fmincon';
options.MS.localOptimizerOptions = optimset('algorithm','interior-point',...
    'display','iter',...
    'GradObj','on',...
    'MaxIter',1000,...
    'TolFun',1e-12,...
    'TolX',1e-10,...
    'MaxFunEvals',10000);
options.MS.obj_type = 'negative log-posterior';

load('erbb_signaling_pnom.mat','pnom');
theta = log10(pnom);
n_par = length(theta);
min_par = theta - 2;
max_par = theta + 3;

par0 = bsxfun(@plus,min_par,bsxfun(@times,max_par - min_par, lhsdesign(n_starts,n_par,'smooth','off')'));
 
switch approach
    case 'standard'
        % nothing to be done
    case 'hierarchical-adjoint'
        n_par = n_par - 3;
        sc.exp_groups.bc_idxs = {1,2,3};
        sc.exp_groups.noise_idxs = {1,2,3};
        
        sc.obs_groups.bc_idcs = {1,2,3};
        sc.obs_groups.b_mode = {'absolute','absolute','absolute'};
        sc.obs_groups.c_mode = {'single','single','single'};
        sc.obs_groups.noise_idxs = {1,2,3};
        sc.obs_groups.noise_mode = {'absolute','absolute','absolute'};
       
        sc.distribution = 'normal';
        
        options.sc = sc;
end

parameters.number = n_par;
parameters.min = min_par;
parameters.max = max_par;
parameters.guess = par0(1:n_par,1:n_starts);

end

