function [] = run_estimation_Zheng(varargin)

approach = varargin{1};
distribution = varargin{2};
noc = varargin{3};

if nargin == 4
    priorsigma2 = varargin{4};
end

addpath(genpath('/home/icb/carolin.loos/PhD/PESTOGit'))
addpath(genpath('/home/icb/carolin.loos/PhD/AMICIGit'))
addpath(genpath('/home/icb/carolin.loos/PhD/HierarchicalOptimization'))

load data_Zheng
[parameters,options] = getParameterOptions_Zheng(approach,noc);

if nargin == 4
    options.llh.prior.flag = true;
    options.llh.prior.mean = -1;
    options.llh.prior.sigma2 = priorsigma2;
else
    options.llh.prior.flag = false;
end

options.llh.distribution = distribution;
options.llh.approach = approach;
options.llh.save_analytical = false;

simfct = @(t,xi,options) simulation_merged_histones_Zheng(t,xi,options);
if options.llh.prior.flag
    if noc
        options.MS.foldername = ['results_Zheng_' approach '_' distribution '_noc_prior_' num2str(priorsigma2)];
    else
        options.MS.foldername = ['results_Zheng_' approach '_' distribution '_prior_' num2str(priorsigma2)];
    end
else
    if noc
        options.MS.foldername = ['results_Zheng_' approach '_' distribution '_noc'];
    else
        options.MS.foldername = ['results_Zheng_' approach '_' distribution];
    end
    
end
if noc
    parameters = getMultiStarts(parameters,@(xi) ...
        logLikelihood_Zheng_noc(xi,D,simfct,options),options.MS);
else
    parameters = getMultiStarts(parameters,@(xi) ...
        logLikelihood_Zheng(xi,D,simfct,options),options.MS);
end
save(options.MS.foldername)

end

%%
% if noc
%     [g,g_fd_f,g_fd_b,g_fd_c] = testGradient(xi,@(xi) logLikelihood_Zheng_noc(xi,D,simfct,options),1e-5)
% else
%     [g,g_fd_f,g_fd_b,g_fd_c] = testGradient(xi,@(xi) logLikelihood_Zheng(xi,D,simfct,options),1e-5)
% end
% [g,g_fd_f,g_fd_b,g_fd_c]
% [g,g_fd_c]