clear all
close all
clc

addpath(genpath('/home/icb/carolin.loos/PhD/PESTOGit'))
addpath(genpath('/home/icb/carolin.loos/PhD/AMICIGit'))
addpath(genpath('/home/icb/carolin.loos/PhD/HierarchicalOptimization'))
addpath(genpath(pwd))

%%
%rmdir('/Users/carolinloos/PhD/Software/AMICIGit/models/Bachmann_JAKSTAT_red','s')
%amiwrap('Bachmann_JAKSTAT_red','Bachmann_JAKSTAT_red_syms',pwd)%

%%
options.llh.approach = 'standard';
options.llh.original = false;
parameters = loadBachmannParameters(options.llh.approach,'reduced_woinit');
parameters.number = numel(parameters.name);
parameters.min = -3*ones(parameters.number,1);
parameters.max = 3*ones(parameters.number,1);
parameters.max(1) = 4; %CISEqc
parameters.max(3) = 12; %CISInh
parameters.max(7) = 4; %EpoRActJAK2
parameters.max(8) = 6; %EpoRCISInh
parameters.max(10) = 9; %JAK2ActEpo
parameters.max(11) = 4; %JAK2EpoRDeaSHP1
parameters.max(20) = 4;
parameters.min(26:56) = -5; %offsets smaller boundary because of changing y = offset + scaling*x to 
                            % y = scaling*(offset + x)
                                                        
xi = loadBestParameter_Bachmann('reduced_woinit'); % reduced = fixed SOCS3RNAEqc and CISRNAEqc and Bachmann_JAKSTAT_red_syms
load data_Bachmann
D = getOffsetScalingStd_Bachmann(D,'reduced_woinit');
% transform data from log scale -> lin scale 
for cond = 1:numel(D)
    D(cond).my(:,[1:10,12:end],:) = 10.^D(cond).my(:,[1:10,12:end],:);
end
D(3).my = D(3).my - 1; %instead of having observable 1 + scaling*...

options.llh.distribution = 'log10-normal';
options.llh.save_analytical = false;
%%
options.llh.reduced_woinit = true;
options.llh.lsqnonlin = 0;
% [ll] = logLikelihood_Bachmann(xi,D,options)
[ll,dll] = neglogLikelihood_Bachmann(xi,D,options);

[g,g_fd_f,g_fd_b,g_fd_c] = testGradient(xi,@(xi) neglogLikelihood_Bachmann(xi,D,options),1e-5);
[g,g_fd_f,g_fd_b,g_fd_c]
%%
options.MS = PestoOptions();
options.MS.localOptimizerOptions = optimset('algorithm','interior-point',...
    'display','iter',...
    'GradObj','on',...
    'MaxIter',5000,...
    'TolFun',1e-10,...
    'TolX',1e-10,...
    'MaxFunEvals',40000,...
    'PrecondBandWidth', inf);

options.MS.comp_type = 'sequential';
options.MS.n_starts = 200;
options.MS.save = true;
options.MS.foldername = 'standard_Bachmann_log10';
options.MS.obj_type = 'negative log-posterior';

load parameter_guesses_Bachmann par0
parameters.guess = par0(1:parameters.number,:);
parameters = getMultiStarts(parameters,@(xi) neglogLikelihood_Bachmann(xi,D,options),options.MS);
save ./results/results_standard_Bachmann_log10

