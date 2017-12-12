addpath(genpath('../'));

clear;

load('data_RafMekErk.mat');
kappa = [[0,0];[0,30];[5,0]];
theta = 0.5*ones(12,1);

rng(0);

n_e = size(D,2); %number of experiments
n_y = size(D(1).my,2); %number of observables
n_r = size(D(1).my,3); %number of replicates

amioptions = amioption();
amioptions.sensi = 2;
amioptions.sensi_meth = 2;

sim = struct([]);
tic
for j=1:n_e
    sim(j).llh=0;
    sim(j).sllh=0;
    sim(j).s2llh=0;
    for k=1:n_r
        clear data_jk;
        tout_j = D(j).t;
        kappa_jk =  [ kappa(j,:) k k ]'; % default values for scalings
        data_jk.Y = D(j).my(:,:,k); % data
        data_jk.Sigma_Y = k*ones(size(data_jk.Y));
        data_jk=amidata(data_jk);
        sol = simulate_RafMekErk_hierarchical_adjoint_reps(tout_j,theta,kappa_jk,data_jk,amioptions);
        if sol.status < 0
            error(['failed to integrate ODE for experiment ' num2str(j)])
        end
        sim(j).llh = sim(j).llh - sol.llh;
        sim(j).sllh = sim(j).sllh - sol.sllh;
        sim(j).s2llh = sim(j).s2llh - sol.s2llh;
    end
end
toc

% and now using amici's ability to handle replicates
sim2 = struct([]);
tic
for j=1:n_e
    clear data_j;
    tout_j = D(j).t;
    kappa_j = zeros(4,n_r);
    for k = 1:n_r
        kappa_j(:,k) = [ kappa(j,:) k k ]'; % default values for scalings
    end
    data_j.Y = D(j).my(:,:,:); % data
    data_j.Sigma_Y = ones(size(data_j.Y));
    for k=1:4
        data_j.Sigma_Y(:,:,k) = k*data_j.Sigma_Y(:,:,k);
    end
    data_j=amidata(data_j);
    sol = simulate_RafMekErk_hierarchical_adjoint_reps(tout_j,theta,kappa_j,data_j,amioptions);
    if sol.status < 0
        error(['failed to integrate ODE for experiment ' num2str(j)])
    end
    sim2(j).llh = - sol.llh;
    sim2(j).sllh = - sol.sllh;
    sim2(j).s2llh = -sol.s2llh;
end
toc

fprintf('Difference in computation:\n');
for j=1:n_e
   fprintf([num2str(sim(j).llh-sim2(j).llh) '  ' mat2str(sim(j).sllh-sim2(j).sllh) '\n']); 
end