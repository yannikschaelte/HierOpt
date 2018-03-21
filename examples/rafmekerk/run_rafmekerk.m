function [  ] = run_rafmekerk(approach)

rng(0);

exdir=fileparts(which('run_rafmekerk'));
addpath(fullfile(exdir,'data'));
addpath(fullfile(exdir,'models'));

[parameters,options] = get_parameters_and_options_rafmekerk(approach);

switch approach
    case 'standard'
        load('data_rafmekerk_noreps.mat','D');
        nllh = @(x) nllh_rafmekerk_standard(x,D);
    case 'hierarchical-forward'
        load('data_rafmekerk.mat','D');
        nllh = @(x) nllh_rafmekerk_hierarchical_forward(x,D,options.sc);
    case 'hierarchical-adjoint'
        load('data_rafmekerk.mat','D');
        nllh = @(x) nllh_rafmekerk_hierarchical_adjoint(x,D,options.sc);
    case 'hierarchical-noreps-forward'
        load('data_rafmekerk_noreps.mat','D');
        nllh = @(x) nllh_rafmekerk_hierarchical_noreps_forward(x,D,options.sc);
    case 'hierarchical-noreps-adjoint'
        load('data_rafmekerk_noreps.mat','D');
        nllh = @(x) nllh_rafmekerk_hierarchical_noreps_adjoint(x,D,options.sc);
    otherwise
        error('approach not recognized');
end

parameters_res = getMultiStarts(parameters,nllh,options.MS);

save(fullfile(exdir,['results_' approach '.mat']));
end

