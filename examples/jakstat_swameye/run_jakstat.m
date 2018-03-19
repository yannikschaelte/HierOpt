function [  ] = run_jakstat(approach)

rng(0);

exdir = fileparts(which('run_jakstat.m'));
addpath(fullfile(exdir,'data'));
addpath(fullfile(exdir,'models'));

load('data_jakstat','D');

[parameters,options] = get_parameters_and_options_jakstat(approach);

switch approach
    case 'standard'
        nllh = @(x) nllh_jakstat_standard(x,D);
    case 'hierarchical-forward'
        nllh = @(x) nllh_jakstat_hierarchical_forward(x,D,options.sc);
    case 'hierarchical-adjoint'
        nllh = @(x) nllh_jakstat_hierarchical_adjoint(x,D,options.sc);
    case 'hierarchical-forward-offsets'
        nllh = @(x) nllh_jakstat_hierarchical_forward_offsets(x,D,options.sc);
    case 'hierarchical-adjoint-offsets'
        nllh = @(x) nllh_jakstat_hierarchical_adjoint_offsets(x,D,options.sc);
    otherwise
        error('Approach not recognized.');
end

parameters_res = getMultiStarts(parameters,nllh,options.MS);

save(fullfile(exdir,['results_' approach '.mat']));

end

