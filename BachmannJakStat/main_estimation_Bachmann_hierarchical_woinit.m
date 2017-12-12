clear all
close all
clc

addpath(genpath('/home/icb/carolin.loos/PhD/PESTOGit'))
addpath(genpath('/home/icb/carolin.loos/PhD/AMICIGit'))
addpath(genpath('/home/icb/carolin.loos/PhD/HierarchicalOptimization'))
addpath(genpath(pwd))
%%
options.llh.approach = 'hierarchical';
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
parameters.min(26:end) = -8; %offsets

xi = loadBestParameter_Bachmann('reduced_woinit');
xi = xi(1:parameters.number)';
xi([5,8,20]) = xi([5,8,20])-0.1;
xi(12) = -3+0.1;
load data_Bachmann
D = getOffsetScalingStd_Bachmann(D);

for cond = 1:numel(D)
    D(cond).my(:,[1:10,12:end],:) = 10.^D(cond).my(:,[1:10,12:end],:);
end
D(3).my = D(3).my - 1;

options.llh.distribution = 'log10-normal';
options.llh.save_analytical = false;
options.llh.reduced_woinit = true;

for i = [1:6,12:20]
    options.llh.obs(i).proportionality = 'single';
    options.llh.obs(i).variance = 'single';
end
for i = 7:11
    options.llh.obs(i).proportionality = 'absolute';
    options.llh.obs(i).variance = 'single';
end
for i = 1:20
    options.llh.obs_groups.proportionality{i} = i;
end
options.llh.obs_groups.variance = {[1,2],[3,19,20],4,[5,6],7,8,9,10,11,[12,13,14,15,16,17],18};
options.llh.exp_groups.proportionality = {1,2,3,[4,5],6,[7,8],[9,10],[11,12],[13,14],...
    [15:19],[20:25],[26:31],[32:36]};

tic
if ~checkValidityOptions(D,options.llh)
    error('options not valid')
end
toc
%% test gradient
options.llh.lsqnonlin = 0;
%xi = getParameterGuesses(parameters,@(xi) logLikelihood_Bachmann(xi,D,options),...
 %    1, parameters.min, parameters.max);
[g,g_fd_f,g_fd_b,g_fd_c] = testGradient(xi,@(xi) neglogLikelihood_Bachmann(xi,D,options),1e-5);
[g,g_fd_f,g_fd_b,g_fd_c]
%%
% options.MS.localOptimizer = 'lsqnonlin';
% options.MS.localOptimizerOptions = optimoptions('lsqnonlin', ...
%     'Algorithm', 'trust-region-reflective', ...
%     'Jacobian', 'on', ...
%     'Display', 'iter-detailed',...
%     'TolFun',1e-16,...
%      'TolX',1e-10,...
%     'MaxIter', 3000);

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
options.MS.foldername = 'hierarchical_Bachmann_log10';
options.MS.obj_type = 'negative log-posterior';

load parameter_guesses_Bachmann par0
parameters.guess = par0(1:parameters.number,:);
parameters = getMultiStarts(parameters,@(xi) neglogLikelihood_Bachmann(xi,D,options),options.MS);
save ./results/results_hierarchical_Bachmann_log10

