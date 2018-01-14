close all
%% VISUALIZATION
%% Font  
setLay.tit.size = 8;
%setLay.tit.font = 'Verdana';
setLay.tit.weight = 'normal';

setLay.lab.size = 8;
%setLay.lab.font = 'Verdana';
setLay.lab.weight = 'normal';

setLay.ax.size = 8;
%setLay.ax.font = 'Verdana';
setLay.lab.weight = 'normal';
title_vec = {'control','30 muM UO126','5 muM Sorafenib'};
%xlab = {'','','','time [h]','time [h]','time [h]'};
xlab = {'','','','','time [h]',''};
ylab = {'relative pMek [UI]','','','relative pErk [UI]','',''};

tsim = (0:0.01:10)';
n_thetaDyn = 12;
kappa = [zeros(1,2);[0,30];[5,0]];
n_kappa = size(kappa,1);
options_simu.sensi = 0;

%%Colors
col_h = [0.2081,0.1663,0.5292]; %color for the hierarchical approach
col_s = [0.9290,0.6940,0.1250]; %color for the standard approach
colors_sh = [col_s;col_h];
color_vp = {[0.9290,0.6940,0.1250], [0.2081,0.1663,0.5292]};

%% Load results of the standard approach
load('results_RafMekErk_standard_normal.mat')
D_s = load('data_RafMekErk_standard');
%% Load results of the hierarchical approach
load('results_RafMekErk_hierarchical_normal.mat')
D_h = load('data_RafMekErk');

[J_opt, grad_opt] = neglogLikelihood_RafMekErk_hierarchical(parameters_hierarchical.MS.par(:,1),kappa,D_h.D,'laplace',1);
load('analytical_results.mat')


%% 'Preprocessing' of the results
for i = 1:n_kappa
    %standard approach
    sol = simulate_RafMekErk(tsim,parameters_standard.MS.par(:,1),kappa(i,:),[],options_simu);    
    bestFit_s(:,:,i) = sol.y;
    
    maxMek_s(i) = max(max(bestFit_s(:,1,i)/bsxfun(@power,10,parameters_standard.MS.par(n_thetaDyn+(1),1))));
    maxErk_s(i) = max(max(bestFit_s(:,2,i)/bsxfun(@power,10,parameters_standard.MS.par(n_thetaDyn+(2),1)))); 
    
    %hierarchical approach
    sol = simulate_RafMekErk_hierarchical(tsim,parameters_hierarchical.MS.par(:,1),kappa(i,:),[],options_simu);
    
    bestFit_h(:,:,i) = sol.y;
    
    maxMek_h(i) = max(max(bestFit_h(:,1,i)));
    maxErk_h(i) = max(max(bestFit_h(:,2,i)));
end

maxMekErk_s = repmat([max(maxMek_s) max(maxErk_s)],1,2);  
maxMekErk_h = repmat([max(maxMek_h) max(maxErk_h)],1,2); 

D_s.D.measurement{1}(:,[3,4,7,8]) = NaN;

%% Best Fit
%figure(1)
for i = 1:n_kappa
    v = figure(i);
    plot(tsim,bestFit_h(:,1,i)/maxMekErk_h(1),'Color',col_h,'LineWidth',1);hold on
    
    if i == 1
        plot(D_h.D(i).t,D_h.D(i).my(:,1,1)/(c(:,1,1)*maxMekErk_h(1)),'og','Markersize',3);hold on
        plot(D_h.D(i).t,D_h.D(i).my(:,1,3)/(c(:,1,3)*maxMekErk_h(1)),'+g','Markersize',3);hold on
    else 
        plot(D_h.D(i).t,D_h.D(i).my(:,1,2)/(c(:,1,2)*maxMekErk_h(1)),'>g','Markersize',3);hold on
        plot(D_h.D(i).t,D_h.D(i).my(:,1,4)/(c(:,1,4)*maxMekErk_h(1)),'xg','Markersize',3);hold on
    end
    
    plot(tsim,bestFit_s(:,1,i)/(bsxfun(@power,10,parameters_standard.MS.par(13,1))*maxMekErk_s(1)),'--','Color',col_s,'LineWidth',1);hold on
  
    if i == 1
        plot(D_s.D.t{i},D_s.D.measurement{i}(:,1)/(bsxfun(@power,10,parameters_standard.MS.par(13,1))*maxMekErk_s(1)),'ok','Markersize',3);hold on
        plot(D_s.D.t{i},D_s.D.measurement{i}(:,5)/(bsxfun(@power,10,parameters_standard.MS.par(17,1))*maxMekErk_s(1)),'+k','Markersize',3);hold on
    else
        plot(D_s.D.t{i},D_s.D.measurement{i}(:,3)/(bsxfun(@power,10,parameters_standard.MS.par(15,1))*maxMekErk_s(1)),'>k','Markersize',3);hold on
        plot(D_s.D.t{i},D_s.D.measurement{i}(:,7)/(bsxfun(@power,10,parameters_standard.MS.par(19,1))*maxMekErk_s(1)),'xk','Markersize',3);hold on
    end
    ylim([-0.2 1.5])
    xlim([-0.1 10.1])
    
     axes = gca;
    set(axes,'Xtick',[0,2,4,6,8,10],'Ytick',[0,0.5,1,1.5],'FontSize',8);
    x0=1;
    y0=1;
    width=85;
    height=70;
    set(gcf,'units','points','position',[x0,y0,width,height]);
    w = figure(i+3);    
    plot(tsim,bestFit_h(:,2,i)/maxMekErk_h(2),'Color',col_h,'LineWidth',1);hold on
    
    if i == 1
        plot(D_h.D(i).t,D_h.D(i).my(:,2,1)/(c(:,2,1)*maxMekErk_h(2)),'og','Markersize',3);hold on
        plot(D_h.D(i).t,D_h.D(i).my(:,2,3)/(c(:,2,3)*maxMekErk_h(2)),'+g','Markersize',3);hold on
    else
        plot(D_h.D(i).t,D_h.D(i).my(:,2,2)/(c(:,2,2)*maxMekErk_h(2)),'>g','Markersize',3);hold on
        plot(D_h.D(i).t,D_h.D(i).my(:,2,4)/(c(:,2,4)*maxMekErk_h(2)),'xg','Markersize',3);hold on
    end
    
    plot(tsim,bestFit_s(:,2,i)/(bsxfun(@power,10,parameters_standard.MS.par(14,1))*maxMekErk_s(2)),'--','Color',col_s,'LineWidth',1);
    
    if i == 1
        plot(D_s.D.t{i},D_s.D.measurement{i}(:,2)/(bsxfun(@power,10,parameters_standard.MS.par(14,1))*maxMekErk_s(2)),'ok','Markersize',3);hold on
        plot(D_s.D.t{i},D_s.D.measurement{i}(:,6)/(bsxfun(@power,10,parameters_standard.MS.par(18,1))*maxMekErk_s(2)),'+k','Markersize',3);hold on
    else
        plot(D_s.D.t{i},D_s.D.measurement{i}(:,4)/(bsxfun(@power,10,parameters_standard.MS.par(16,1))*maxMekErk_s(2)),'>k','Markersize',3);hold on
        plot(D_s.D.t{i},D_s.D.measurement{i}(:,8)/(bsxfun(@power,10,parameters_standard.MS.par(20,1))*maxMekErk_s(2)),'xk','Markersize',3);hold on
    end
    ylim([-0.2 1.5])
    xlim([-0.1 10.1])
    
    axes = gca;
    set(axes,'Xtick',[0,2,4,6,8,10],'Ytick',[0,0.5,1,1.5],'FontSize',8);
    x0=1;
    y0=1;
    width=80;
    height=70;
    set(gcf,'units','points','position',[x0,y0,width,height]);
end

%% Convergence
starts_h = parameters_hierarchical.MS.logPost;
starts_h(abs(starts_h-starts_h(1))>1e-3) = NaN;
starts_h(~isnan(starts_h)) = 1;
n_conv_starts_s = nansum(starts_h);
c = figure(7);
plot(-parameters_hierarchical.MS.logPost,'o-','Color',col_h,'MarkerSize',1);
hold on
plot(-parameters_standard.MS.logPost,'o-','Color',col_s,'MarkerSize',1);
hold on
plot((1:2:n_conv_starts_s),-parameters_hierarchical.MS.logPost(1:2:n_conv_starts_s),'o-','Color',col_h,'MarkerSize',1);
xlim([0,250]);
ylim([45,70]);
set(gca,'FontSize',setLay.ax.size)
x0=1;
y0=1;
width=250;
height=100;
set(gcf,'units','points','position',[x0,y0,width,height]);

%% Boxplots cpu time
approach = {'standard','hierarchical'};
cpu_time_n = parameters_standard.MS.t_cpu;
cpu_time_n(cpu_time_n==0) = NaN;
cpu_time_a = parameters_hierarchical.MS.t_cpu;
cpu_time_a(cpu_time_a==0) = NaN;
bp = figure(8);
boxplot([log10(cpu_time_n),log10(cpu_time_a)],approach,'OutlierSize',4,'Symbol','k.');hold on;
boxes = findobj(gca,'Tag','Box');
median = findobj(gca,'Tag','Median');
patch(get(boxes(1),'XData'),get(boxes(1),'YData'),col_h,'FaceAlpha',1); hold on;
patch(get(boxes(2),'XData'),get(boxes(2),'YData'),col_s,'FaceAlpha',1); hold on;
plot(get(median(1),'XData'),get(median(1),'YData'),'k-','Linewidth',0.5); hold on;
plot(get(median(2),'XData'),get(median(2),'YData'),'k-','Linewidth',0.5);
ylim([0,3]);
set(gca,'YTick',(0:1:3));
set(gca,'XTickLabel',{'standard','hierarchical'},'YTickLabel',{'10^0' '10^1' '10^2' '10^3'},'TickLabelInterpreter', 'tex');
set(gca,'FontSize',setLay.ax.size)
x0=1;
y0=1;
width = 140;
height = 100;
set(gcf,'units','points','position',[x0,y0,width,height]);

%% Cpu time/converged start
starts_n = parameters_standard.MS.logPost;
starts_n(abs(starts_n-starts_n(1))>10e-3) = NaN;
starts_n(~isnan(starts_n)) = 1;
conv_starts_n = nansum(starts_n);
starts_a = parameters_hierarchical.MS.logPost;
starts_a(abs(starts_a-starts_a(1))>10e-3) = NaN;
starts_a(~isnan(starts_a)) = 1;
conv_starts_a = nansum(starts_a);
cpu_time_n = sum(parameters_standard.MS.t_cpu);
cpu_time_a = sum(parameters_hierarchical.MS.t_cpu);
vs_n = [log10(cpu_time_n/conv_starts_n),0];
vs_a = [0,log10(cpu_time_a/conv_starts_a)];
v = figure(9);
b_n = bar(vs_n,'BarWidth',0.5,'FaceColor',col_s);
hold on
b_a = bar(vs_a,'BarWidth',0.5,'FaceColor',col_h);
set(gca,'YTick',(0:1:3));
set(gca, 'XTickLabel',{'standard','hierarchical'},'YTickLabel',{'10^0' '10^1' '10^2' '10^3'},'TickLabelInterpreter', 'tex');
set(gca,'FontSize',setLay.ax.size);
 x0=1;
y0=1;
width = 140;
height = 100;
set(gcf,'units','points','position',[x0,y0,width,height]);

%%
load('results_RafMekErk_standard_normal.mat')
load('results_RafMekErk_hierarchical_normal.mat')
cpu_time_n = parameters_standard.MS.t_cpu;
cpu_time_n(cpu_time_n==0) = NaN;
cpu_time_a = parameters_hierarchical.MS.t_cpu;
cpu_time_a(cpu_time_a==0) = NaN;

cpu_time_n = sum(parameters_standard.MS.t_cpu);%
cpu_time_a = sum(parameters_hierarchical.MS.t_cpu);%
starts_n = parameters_standard.MS.logPost;
starts_n(abs(starts_n-starts_n(1))>1e-3) = NaN;
starts_n(~isnan(starts_n)) = 1;
%conv_starts_n = nansum(starts_n)
conv_starts_n = sum(2*(parameters_standard.MS.logPost-parameters_standard.MS.logPost(1))>-icdf('chi2',0.95,1));
conv_starts_n/sum(~isnan(parameters_standard.MS.logPost))
cpu_time_n/conv_starts_n

starts_a = parameters_hierarchical.MS.logPost;
starts_a(abs(starts_a-starts_a(1))>1e-3) = NaN;
starts_a(~isnan(starts_a)) = 1;
%conv_starts_a = nansum(starts_a);
conv_starts_a = sum(2*(parameters_hierarchical.MS.logPost-parameters_hierarchical.MS.logPost(1))>-icdf('chi2',0.95,1));

conv_starts_a/sum(~isnan(parameters_hierarchical.MS.logPost))
cpu_time_a/conv_starts_a
