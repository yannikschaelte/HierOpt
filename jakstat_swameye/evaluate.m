close all;

load results_standard.mat
pam_std = parameters_res;
load results_hierarchical.mat
pam_h = parameters_res;

% simulate
% load data/data_RafMekErk_standard.mat
% i = 1;
% u = D.conditions;
% t_m = D.t{i};
% n_t = length(t_m);
% options_simu.sensi = 0;
% sol = simulate_RafMekErk(t_m,pam_std.MS.par,u(i,:),[],options_simu);
% figure;
% plot(sol.y(1,:));

load data/data_jakstat.mat

arr_pams = {pam_std,pam_h};%,pam_ha,pam_ha2};
identifiers = {'standard','hierarchical'};
for j = 1:length(arr_pams)
    figure;
    fig = plotSimple(arr_pams{j});
    title(['All estimates (' identifiers{j} ')']);
%     plotMultiStartDiagnosis(arr_pams{j});
    t(j) = median(arr_pams{j}.MS.t_cpu);
    fprintf('%d: %.15f\n',j,t(j));
    arr_mean(j) = mean(arr_pams{j}.MS.t_cpu);
    arr_std(j) = std(arr_pams{j}.MS.t_cpu);
    t_cpu(:,j) = arr_pams{j}.MS.t_cpu;
    nConvergedStarts(1,j) = sum(arr_pams{j}.MS.logPost > 73);
end

figure;
b = bar(nConvergedStarts/100*100);
xticklabels({'standard','hierarchical'});
ylabel('converged starts [%]');
title('Converged starts');

figure;
boxplot(t_cpu,'DataLim',[0,50]);
xticklabels({'standard','hierarchical'});
ylabel('time [s]');
title('CPU time per start');

figure;
plot(t,'o');
ylim([0,Inf]);