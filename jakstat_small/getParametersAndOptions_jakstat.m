function [parameters,options] = getParametersAndOptions_jakstat(approach)

nStarts = 50;

options.MS = PestoOptions();
options.MS.n_starts = nStarts; %actually 500
options.MS.mode = 'text';
options.MS.localOptimizer = 'fmincon';
options.MS.localOptimizerOptions = optimset('algorithm','interior-point',...
    'display','iter',...
    'GradObj','on',...
    'MaxIter',2000,...
    'TolFun',1e-10,...
    'TolX',1e-10,...
    'MaxFunEvals',2000,...
    'PrecondBandWidth', inf);
options.MS.obj_type = 'negative log-posterior';

nPar = 17;
minPar = -5*ones(nPar,1);
maxPar = 3*ones(nPar,1);
maxPar(4)  =  6;
maxPar(2)  =  6;
minPar(10) = -6;
minPar(4)  = -3;
minPar(2)  = -3;

rng(1);
par0 = bsxfun(@plus,minPar,bsxfun(@times,maxPar - minPar, lhsdesign(nStarts,nPar,'smooth','off')'));
   
switch approach
    case 'standard'
        nPar = 17;
        
    case 'hierarchical-adjoint'
        nPar = 12;
        options.sc.obs(1).variance = 'multiple';
        options.sc.obs(1).proportionality = 'multiple';
        options.sc.obs(2).variance = 'multiple';
        options.sc.obs(2).proportionality = 'multiple';
        options.sc.obs(3).variance = 'multiple';
        options.sc.obs(3).proportionality = 'absolute';
        options.sc.obs_groups.variance = {1,2,3};
        options.sc.obs_groups.proportionality = {1,2,3};
    
    case 'hierarchical'
        nPar = 12;
        options.sc.obs(1).variance = 'multiple';
        options.sc.obs(1).proportionality = 'multiple';
        options.sc.obs(2).variance = 'multiple';
        options.sc.obs(2).proportionality = 'multiple';
        options.sc.obs(3).variance = 'multiple';
        options.sc.obs(3).proportionality = 'absolute';
        options.sc.obs_groups.variance = {1,2,3};
        options.sc.obs_groups.proportionality = {1,2,3};
        
    case 'hierarchical-offsets'
        nPar = 10;
        options.sc.obs(1).variance = 'multiple';
        options.sc.obs(1).proportionality = 'multiple';
        options.sc.obs(2).variance = 'multiple';
        options.sc.obs(2).proportionality = 'multiple';
        options.sc.obs(3).variance = 'multiple';
        options.sc.obs(3).proportionality = 'absolute';
        options.sc.obs_groups.variance = {1,2,3};
        options.sc.obs_groups.proportionality = {1,2,3};
        
    case 'hierarchical-adjoint-offsets'
        nPar = 10;
        options.sc.obs(1).variance = 'multiple';
        options.sc.obs(1).proportionality = 'multiple';
        options.sc.obs(2).variance = 'multiple';
        options.sc.obs(2).proportionality = 'multiple';
        options.sc.obs(3).variance = 'multiple';
        options.sc.obs(3).proportionality = 'absolute';
        options.sc.obs_groups.variance = {1,2,3};
        options.sc.obs_groups.proportionality = {1,2,3};
end

parameters.number = nPar;
parameters.min = minPar(1:nPar,1);
parameters.max = maxPar(1:nPar,1);
parameters.guess = par0(1:nPar,1:nStarts);

end

