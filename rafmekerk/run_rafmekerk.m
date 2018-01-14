function [  ] = run_rafmekerk(approach)

kappa = [zeros(1,2);[0,30];[5,0]];

exdir=fileparts(which('run_rafmekerk'));
[parameters,options] = getParametersAndOptions_rafmekerk(approach);

switch approach
    case 'standard'
        load data_RafMekErk_standard.mat
        nllh = @(x) nllh_rafmekerk_standard(x,kappa,D);
    case 'hierarchical'
        load data_RafMekErk.mat
        nllh = @(x) nllh_rafmekerk_hierarchical(x,kappa,D,options.sc);
    case 'hierarchical-adjoint'
        load data_RafMekErk.mat
        nllh = @(x) nllh_rafmekerk_hierarchical_adjoint(x,kappa,D,options.sc);
    otherwise
        error('approach not recognized');
end

parameters_res = getMultiStarts(parameters,nllh,options.MS);

save(fullfile(exdir,['results_' approach '.mat']));
end

