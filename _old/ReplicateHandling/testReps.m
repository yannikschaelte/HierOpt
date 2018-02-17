addpath(genpath('..'));

rng(0);

maxr = 100;
arr_r = 1:maxr;
arr_time1 = zeros(maxr,1);
arr_time2 = zeros(maxr,1);
%%

for j=1:maxr
    load D_jakstat_pesto.mat
    
    [n_t, n_y] = size(D);
    n_r = arr_r(j);
    
    datatable = zeros(n_t,n_y,n_r);
    
    for it = 1:n_t
        for iy = 1:n_y
            for ir = 1:n_r
                datatable(it,iy,ir) = D(it,iy) + 1e-2*randn(1);
            end
        end
    end
    
    clear D;
    
    D.my = datatable;
    D.kappa = [1.4,0.45];
    D.t = 1:n_t;
    
    theta = rand(17,1);
    
    %%
    
    amioptions = amioption();
    amioptions.sensi = 1;
    amioptions.sensi_meth = 2;
    
    sim = struct();
    sim.llh=0;
    sim.sllh=0;
    sim.s2llh=0;
    starttime = cputime;
    for k=1:n_r
        clear data_jk;
        tout_j = D.t;
        kappa_j = D.kappa; % default values for scalings
        data_jk.Y = D.my(:,:,k); % data
        data_jk.Sigma_Y = nan(size(data_jk.Y));
        data_jk=amidata(data_jk);
        sol = simulate_jakstat_pesto(tout_j,theta,kappa_j,data_jk,amioptions);
        if sol.status < 0
            error('failed to integrate ODE for experiment')
        end
        sim.llh = sim.llh - sol.llh;
        sim.sllh = sim.sllh - sol.sllh;
        sim.s2llh = sim.s2llh - sol.s2llh;
    end
    arr_time1(j) = cputime - starttime;
    
    % and now using amici's ability to handle replicates
    sim2 = struct();
    starttime = cputime;
    clear data_j;
    tout_j = D.t;
    kappa_j = D.kappa; % default values for scalings
    data_j.Y = D.my(:,:,:); % data
    data_j.Sigma_Y = nan(size(data_j.Y));
    data_j=amidata(data_j);
    sol = simulate_jakstat_pesto(tout_j,theta,kappa_j,data_j,amioptions);
    if sol.status < 0
        error('failed to integrate ODE for experiment')
    end
    sim2.llh = - sol.llh;
    sim2.sllh = - sol.sllh;
    sim2.s2llh = -sol.s2llh;
    arr_time2(j) = cputime - starttime;
    
end

figure;
hold on;
plot(arr_r,arr_time1,'-','DisplayName','old');
plot(arr_r,arr_time2,'-','DisplayName','new');
xlabel('replicates');
ylabel('time [s]');
legend('show','Location','northeastoutside');
title('Efficient adjoints with replicates');
figure;
plot(arr_r,arr_time1./arr_time2);