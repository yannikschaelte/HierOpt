function [parameters,options] = getParameterOptions_JakStat(approach)

options.MS = PestoOptions();
options.MS.n_starts = 20; %actually 500
options.MS.mode = 'text';
options.MS.localOptimizer = 'fmincon';
options.MS.localOptimizerOptions = optimset('algorithm','interior-point',...
    'display','iter',...
    'GradObj','on',...
    'MaxIter',5000,...
    'TolFun',1e-10,...
    'TolX',1e-10,...
    'MaxFunEvals',40000,...
    'PrecondBandWidth', inf);
options.MS.obj_type = 'negative log-posterior';


load parameter_guesses_JakStat par0

switch approach
    case {'hierarchical','hierarchical-adjoint','hierarchical-adjoint-offsets'}
        parameters.name = {'log_{10}(p1)','log_{10}(p2)','log_{10}(p3)','log_{10}(p4)',...
            'log_{10}(sp1)','log_{10}(sp2)','log_{10}(sp3)','log_{10}(sp4)','log_{10}(sp5)',...
            'log_{10}(offset_{tSTAT})','log_{10}(offset_{pSTAT})'};
%         parameters.guess = par0(1:length(parameters.name),1:options.MS.n_starts);
        options.llh.obs(1).variance = 'multiple';
        options.llh.obs(1).proportionality = 'multiple';
        options.llh.obs(2).variance = 'multiple';
        options.llh.obs(2).proportionality = 'multiple';
        options.llh.obs(3).variance = 'multiple';
        options.llh.obs(3).proportionality = 'absolute';
        options.llh.obs_groups.variance = {1,2,3};
        options.llh.obs_groups.proportionality = {1,2,3};
        
    case 'standard'
        
        parameters.name = {'log_{10}(p1)','log_{10}(p2)','log_{10}(p3)','log_{10}(p4)',...
            'log_{10}(sp1)','log_{10}(sp2)','log_{10}(sp3)','log_{10}(sp4)','log_{10}(sp5)',...
            'log_{10}(offset_{tSTAT})','log_{10}(offset_{pSTAT})',...
            'log_{10}(scale_{tSTAT})','log_{10}(scale_{pSTAT})',...
            'log_{10}(\sigma_{pSTAT})','log_{10}(\sigma_{tSTAT})','log_{10}(\sigma_{pEpoR})'};
%         parameters.guess = par0(:,1:options.MS.n_starts);
end

parameters.number = length(parameters.name);
parameters.min = -5*ones(parameters.number,1);
parameters.max = 3*ones(parameters.number,1);
parameters.max(4) = 6;
parameters.max(2) = 6;
parameters.min(9) = -6;
parameters.min(4) = -3;
parameters.min(2) = -3;