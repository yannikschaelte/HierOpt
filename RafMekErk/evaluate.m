close all;

load results_RafMekErk_standard_normal_500.mat
pam_std = parameters_res;
load results_RafMekErk_hierarchical_normal_500.mat
pam_h = parameters_res;
load results_RafMekErk_hierarchical-adjoint_normal_500.mat
pam_ha = parameters_res;
load results_RafMekErk_hierarchical-adjoint2_normal_500.mat
pam_ha2 = parameters_res;

arr_pams = {pam_std,pam_h,pam_ha,pam_ha2};
for j = 1:length(arr_pams)
    plotMultiStarts(arr_pams{j});
    plotMultiStartDiagnosis(arr_pams{j});
    t(j) = median(arr_pams{j}.MS.t_cpu);
    fprintf('%d: %.15f\n',j,t(j));
end

figure;
plot(t,'o');
ylim([0,Inf]);