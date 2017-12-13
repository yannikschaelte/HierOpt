function [  ] = run_jakstat(approach)

addpath(genpath('..'));

load('data_jakstat','D');
kappa = [1.4; 0.45]; % omega_cyt, omega_nuc

[parameters,options] = getParametersAndOptions_jakstat(approach);

switch approach
    case 'standard'
        nllh = @(x) nllh_jakstat_standard(x,kappa,D);
    case 'hierarchical'
        nllh = @(x) nllh_jakstat_hierarchical(x,kappa,D,options.sc);
    case 'hierarchical-adjoint'
        nllh = @(x) nllh_jakstat_hierarchical_adjoint(x,kappa,D,options.sc);
end

parameters_res = getMultiStarts(parameters,nllh,options.MS);

save(['results_' approach '.mat']);

end

