function [] = run_estimation_RafMekErk(varargin)

% addpath(genpath('/home/icb/carolin.loos/PhD/PESTOGit'))
% addpath(genpath('/home/icb/carolin.loos/PhD/AMICIGit'))
% addpath(genpath('/home/icb/carolin.loos/PhD/HierarchicalOptimization'))
addpath(genpath('..'));

rng(0);

if nargin == 0
    approach = 'hierarchical';
    distribution = 'normal';
	use_prior = false;
else
    approach = varargin{1};
    distribution = varargin{2};
    if nargin == 3
        priorsigma2 = varargin{3};
        use_prior = true;
	else
		use_prior = false;  
	end
end

switch approach
    case 'standard'
        load('data_RafMekErk_standard')
    case 'hierarchical'
        load('data_RafMekErk')
end

[parameters,options] = getParameterOptions_RafMekErk(approach);

if use_prior
    options.MS.foldername = ['results_RafMekErk_' approach '_' distribution '_prior' num2str(priorsigma2)];
else
    options.MS.foldername = ['results_RafMekErk_' approach '_' distribution];
end

switch approach
    case 'hierarchical'
        load data_RafMekErk.mat
        kappa = [zeros(1,2);[0,30];[5,0]];
        if use_prior
            parameters_res = getMultiStarts(parameters,@(xi) ...
                neglogLikelihood_RafMekErk_hierarchical(xi,kappa,D,distribution,0,priorsigma2),...
                options.MS);
        else
            parameters_res = getMultiStarts(parameters,@(xi) ...
                neglogLikelihood_RafMekErk_hierarchical(xi,kappa,D,distribution,0),...
                options.MS);
        end
    case 'hierarchical-adjoint'
        load data_RafMekErk.mat
        kappa = [[0,0];[0,30];[5,0]];
        parameters_res = getMultiStarts(parameters,@(xi) ...
            neglogLikelihood_RafMekErk_hierarchical_adjoint(xi,kappa,D,distribution,0),...
            options.MS);
    case 'hierarchical-adjoint2'
        load data_RafMekErk_standard.mat
        parameters_res = getMultiStarts(parameters,@(xi) ...
            neglogLikelihood_RafMekErk_hierarchical_adjoint2(xi,D,distribution,0),...
            options.MS);
    case 'hierarchical-adjoint2-offsets'
        load data_RafMekErk_standard.mat
        parameters_res = getMultiStarts(parameters,@(xi) ...
            neglogLikelihood_RafMekErk_hierarchical_adjoint2_offsets(xi,D,distribution,0),...
            options.MS);
    case 'hierarchical-adjoint-reps'
        load data_RafMekErk.mat
        kappa = [[0,0];[0,30];[5,0]];
        parameters_res = getMultiStarts(parameters,@(xi) ...
            neglogLikelihood_RafMekErk_hierarchical_adjoint_reps(xi,kappa,D,distribution,0),...
            options.MS);
    case 'standard'
        load data_RafMekErk_standard.mat
        if use_prior
            parameters_res = getMultiStarts(parameters,@(xi) ...
                neglogLikelihood_RafMekErk_standard(xi,D,distribution,0,priorsigma2),options.MS);
        else
            parameters_res = getMultiStarts(parameters,@(xi) ...
                neglogLikelihood_RafMekErk_standard(xi,D,distribution),options.MS);
        end
end
save(options.MS.foldername)
end

