function [ ] = run_erbb(approach)

exdir = fileparts(which('run_erbb.m'));

D = get_data_erbb();

[parameters,options] = get_parameters_and_options_erbb(approach);

switch approach
    case 'standard'
        nllh = @(x) nllh_erbb_standard(x,D);
end

parameters_res = getMultiStarts(parameters, nllh, options.MS);

save(fullfile(exdir,['results_' approach '.mat']));

end