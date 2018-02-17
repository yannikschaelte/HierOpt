function [] = run_estimation_JakStat(approach,distribution)

% addpath(genpath('/home/icb/carolin.loos/PhD/PESTOGit'))
% addpath(genpath('/home/icb/carolin.loos/PhD/AMICIGit'))
% addpath(genpath('/home/icb/carolin.loos/PhD/HierarchicalOptimization'))

load data_JakStat.mat
[parameters,options] = getParameterOptions_JakStat(approach);

kappa(1) = 1.4;% Omega_cyt
kappa(2) = 0.45;% Omega_nuc
kappa(3) = 1; % init_STAT

options.llh.distribution = distribution;
options.llh.approach = approach;
options.llh.save_analytical = false;

options.MS.foldername = ['results_SmallJakStat_' approach '_' distribution];

switch approach
    case 'hierarchical-adjoint'
        parameters_res = getMultiStarts(parameters,@(xi)...
            neglogLikelihood_JakStat_hierarchical_adjoint(xi,kappa,D,distribution,options),options.MS);
    case 'hierarchical-adjoint-offsets'
        parameters_res = getMultiStarts(parameters,@(xi)...
            neglogLikelihood_JakStat_hierarchical_adjoint_offsets(xi,kappa,D,distribution,options),options.MS);
    case {'standard','hierarchical'}
        parameters_res = getMultiStarts(parameters,@(xi) ...
            neglogLikelihood_JakStat(xi,kappa,D,distribution,options),options.MS);
    otherwise
        error('Approach not recognized.');
end

save(options.MS.foldername);
end