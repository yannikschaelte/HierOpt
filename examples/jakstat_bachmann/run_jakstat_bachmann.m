function [ ] = run_jakstat_bachmann(approach)

if nargin == 0
    approach = 'standard';
end % TODO remove

exdir = fileparts(which('run_jakstat_bachmann'));

[parameters,options] = getParametersAndOptions_jakstat_bachmann(approach);
D = getData();

nllh = @(x) nllh_jakstat_bachmann_standard(x,D);

parameters_res = getMultiStarts(parameters,nllh,options.MS);

save(fullfile(exdir,['results_' approach '.mat']));

end

