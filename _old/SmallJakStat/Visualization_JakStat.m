clear all
close all
clc

%% Font  
setLay.tit.size = 14;
setLay.tit.font = 'Verdana';
setLay.tit.weight = 'normal';

setLay.lab.size = 14;
setLay.lab.font = 'Verdana';
setLay.lab.weight = 'bold';

setLay.ax.size = 14;
setLay.ax.font = 'Verdana';

%% Load data
load('results_JakStat_standard.mat')
parameters_standard = parameters;
load('results_JakStat_hierarchical.mat')
parameters_hierarchical = parameters;
options.llh.save_analytical = true;
neglogLikelihood_JakStat(parameters_hierarchical.MS.par(:,1),kappa,D,distr,options)
load('analytical_results.mat')
load('data_JakStat.mat')

%% Simulation
tsim = (0:0.5:60)';
options.ami.sensi = 0;

sol_n = simulate_JakStat(tsim,parameters_standard.MS.par(1:14,1),kappa,[]);
sol_r = simulate_JakStat_hierarchical(tsim,parameters_hierarchical.MS.par(1:11,1),kappa,[],options.ami);

%% Visualization
%Colors
col_h = [0.2081,0.1663,0.5292]; %color for the hierarchical approach
col_s = [0.9290,0.6940,0.1250]; %color for the standard approach
colors_sh = [col_s;col_h];
color_vp = {[0.9290,0.6940,0.1250], [0.2081,0.1663,0.5292]};

%% Fit
y1 = figure(4);
plot(tsim,c(1)*sol_r.y(:,1),'Color',col_h,'LineWidth',2);hold on
plot(tsim,sol_n.y(:,1),'--','Color',col_s,'LineWidth',2);hold on
plot(D.t,D.my(:,1),'.','MarkerSize',15,'MarkerEdgeColor','black','MarkerFaceColor','black');
xlim([0 60]);
ylim([0 1.1]);
set(gca,'FontSize',setLay.ax.size)
set(gca,'XTick',(0:10:60),'XtickLabel',{''});
pos = get(y1,'Position');
set(y1,'Position',[pos(1) pos(2) 0.5*pos(3) 0.25*pos(4)]); 
y2 = figure(5);
%subplot(3,1,2)
plot(tsim,c(2)*sol_r.y(:,2),'Color',col_h,'LineWidth',2);hold on
plot(tsim,sol_n.y(:,2),'--','Color',col_s,'LineWidth',2);hold on
plot(D.t,D.my(:,2),'.','MarkerSize',15,'MarkerEdgeColor','black','MarkerFaceColor','black');
xlim([0 60]);
ylim([0 1.1]);
set(gca,'FontSize',setLay.ax.size)
set(gca,'XTick',(0:10:60),'XtickLabel',{''});
pos = get(y2,'Position');
set(y2,'Position',[pos(1) pos(2) 0.5*pos(3) 0.25*pos(4)]); 

y3 = figure(6);
%subplot(3,1,3)
plot(tsim,sol_r.y(:,3),'Color',col_h,'LineWidth',2);hold on
plot(tsim,sol_n.y(:,3),'--','Color',col_s,'LineWidth',2);hold on 
plot(D.t,D.my(:,3),'.','MarkerSize',15,'MarkerEdgeColor','black','MarkerFaceColor','black');
xlim([0 60]);
ylim([0 1.1]);
set(gca,'FontSize',setLay.ax.size)
set(gca,'XTick',(0:10:60));
pos = get(y3,'Position');
set(y3,'Position',[pos(1) pos(2) 0.5*pos(3) 0.25*pos(4)]); 

%% Cpu Time
% Boxplots cpu time
approach = {'standard','hierarchical'};
cpu_time_n = parameters_standard.MS.t_cpu;
cpu_time_n(cpu_time_n==0) = NaN;
cpu_time_a = parameters_hierarchical.MS.t_cpu;
cpu_time_a(cpu_time_a==0) = NaN;
bp = figure(1);
%boxplot([log10(cpu_time_n),log10(cpu_time_a)],approach,'Color', colors_sh,'OutlierSize',3,'Symbol','o');
boxplot([log10(cpu_time_n),log10(cpu_time_a)],approach,'OutlierSize',8,'Symbol','k.');hold on;
boxes = findobj(gca,'Tag','Box');
median = findobj(gca,'Tag','Median');
patch(get(boxes(1),'XData'),get(boxes(1),'YData'),col_h,'FaceAlpha',1); hold on;
patch(get(boxes(2),'XData'),get(boxes(2),'YData'),col_s,'FaceAlpha',1); hold on;
plot(get(median(1),'XData'),get(median(1),'YData'),'k-','Linewidth',0.5); hold on;
plot(get(median(2),'XData'),get(median(2),'YData'),'k-','Linewidth',0.5);
ylim([-2.1,4.0]);
set(gca,'YTick',(-2:1:4));
set(gca,'XTickLabel',{'',''},'YTickLabel',{'10^{-2}' '10^{-1}' '10^0' '10^1' '10^2' '10^3' '10^4'},'TickLabelInterpreter', 'tex');
pos = get(bp,'Position');
set(bp,'Position',[pos(1) pos(2) 0.45*pos(3) 0.92*pos(4)]); 

%% Cpu time/converged start
starts_n = parameters_standard.MS.logPost;
starts_n(abs(starts_n-starts_n(1))>1e-3) = NaN;
starts_n(~isnan(starts_n)) = 1;
%conv_starts_n = nansum(starts_n)
conv_starts_n = sum(2*(parameters_standard.MS.logPost-parameters_standard.MS.logPost(1))>-icdf('chi2',0.95,1));

conv_starts_n/sum(~isnan(parameters_standard.MS.logPost))

starts_a = parameters_hierarchical.MS.logPost;
starts_a(abs(starts_a-starts_a(1))>1e-3) = NaN;
starts_a(~isnan(starts_a)) = 1;
%conv_starts_a = nansum(starts_a);
conv_starts_a = sum(2*(parameters_hierarchical.MS.logPost-parameters_hierarchical.MS.logPost(1))>-icdf('chi2',0.95,1));
conv_starts_a/sum(~isnan(parameters_hierarchical.MS.logPost))

cpu_time_n = sum(parameters_standard.MS.t_cpu);%
cpu_time_a = sum(parameters_hierarchical.MS.t_cpu);%
vs_n = [log10(cpu_time_n/conv_starts_n),0];%
vs_a = [0,log10(cpu_time_a/conv_starts_a)];

v = figure(2);
%subplot(2,2,4)
b_n = bar(vs_n,'BarWidth',0.5,'FaceColor',col_s);
hold on
b_a = bar(vs_a,'BarWidth',0.5,'FaceColor',col_h);
%ylim([0,3.5]);
set(gca,'YTick',(0:1:3));
set(gca, 'FontSize', setLay.ax.size);%,'FontName',setLay.lab.font);
set(gca, 'XTickLabel',{'standard','hierarchical'},'YTickLabel',{'10^0' '10^1' '10^2' '10^3'},'TickLabelInterpreter', 'tex');
pos = get(v,'Position');
set(v,'Position',[pos(1) pos(2) 0.45*pos(3) 0.5*pos(4)]);

%% Convergence plot
starts_h = parameters_standard.MS.logPost;
starts_h(abs(starts_h-starts_h(1))>1e-1) = NaN;
starts_h(~isnan(starts_h)) = 1;
n_conv_starts_s = nansum(starts_h);
c = figure(3);
plot(-parameters_hierarchical.MS.logPost,'o-','Color',col_h);
hold on
plot(-parameters_standard.MS.logPost,'o-','Color',col_s);
hold on
plot((1:2:n_conv_starts_s),-parameters_hierarchical.MS.logPost(1:2:n_conv_starts_s),'o-','Color',col_h);
xlim([0,65]);
ylim([-100,130]);
set(gca,'FontSize',setLay.ax.size)%,'FontName',setLay.lab.font)
pos = get(c,'Position');
set(c,'Position',[pos(1) pos(2) 0.5*pos(3) 0.5*pos(4)]);


%%
load('results_JakStat_standard_laplace.mat')
load('results_JakStat_hierarchical_laplace.mat')
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

