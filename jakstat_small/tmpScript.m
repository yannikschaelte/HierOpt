clear;
load results_hierarchical-offsets.mat;
theta = parameters_res.MS.par(:,1);
amiOptions.sensi = 1;
amiOptions.sensi_meth = 'forward';

sim = simulate_jakstat_hierarchical_offsets(D.t,theta,kappa,[],amiOptions);

options.sc.obs(1).variance = 'multiple';
        options.sc.obs(1).proportionality = 'multiple';
        options.sc.obs(2).variance = 'multiple';
        options.sc.obs(2).proportionality = 'multiple';
        options.sc.obs(3).variance = 'multiple';
        options.sc.obs(3).proportionality = 'absolute';
        options.sc.obs_groups.variance = {1,2,3};
        options.sc.obs_groups.proportionality = {1,2,3};
scOptions = options.sc;
        
[c,sigma2,~,~,b,~] = opt_scalings_normal(sim,D,scOptions);
[nllh] = opt_nllh(D,sim,b,c,sigma2)
thetaComplete = [theta;b(2:-1:1)';c(2:-1:1)';sigma2(1,:)'];
D.Y = D.my;
simComplete = simulate_jakstat(D.t,thetaComplete,kappa,D,amiOptions);
-simComplete.llh

load results_hierarchical.mat
theta2 = parameters_res.MS.par(:,1);
sim = simulate_jakstat_hierarchical(D.t,theta2,kappa,[],amiOptions);
[c,sigma2] = scalings(sim,D,scOptions);
[nllh] = nlLH_fgh(sim,D,'normal',scOptions)
