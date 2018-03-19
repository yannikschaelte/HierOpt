function [ ] = run_erbb(approach)

rng(0);

exdir = fileparts(which('run_erbb.m'));
addpath(fullfile(exdir,'data'));
addpath(fullfile(exdir,'models'));

D = get_data_erbb();

[parameters,options] = get_parameters_and_options_erbb(approach);

switch approach
    case 'standard'
        nllh = @(x) nllh_erbb_standard(x,D);
    case 'hierarchical-adjoint'
        nllh = @(x) nllh_erbb_hierarchical_adjoint(x,D,options.sc);
    otherwise
        error('Approach not recognized.');
end

parameters_res = getMultiStarts(parameters, nllh, options.MS);

save(fullfile(exdir,['results_' approach '.mat']));

end