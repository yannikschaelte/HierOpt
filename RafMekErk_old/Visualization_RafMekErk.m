%% VISUALISATION
%% Font  
setLay.tit.size = 14;
%setLay.tit.font = 'Verdana';
setLay.tit.weight = 'normal';

setLay.lab.size = 14;
%setLay.lab.font = 'Verdana';
setLay.lab.weight = 'normal';

setLay.ax.size = 14;
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
% % col_h = [0,0.4470,0.7410]; %color for the hierarchical approach
% % col_s = [0.8500,0.3250,0.0980]; % color for the standard approach
% % colors_sh = [col_s;col_h];
% % color_vp = {[0.8500,0.3250,0.0980], [0,0.4470,0.7410]};
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
    %[status,~,~,y_i] = simulate_RafMekErklin_ana(tsim,parameters_standard.MS.par(:,1),kappa(i,:));
    
    %bestFit_s(:,:,i) = y_i;
    bestFit_s(:,:,i) = sol.y;
    
    maxMek_s(i) = max(max(bestFit_s(:,1,i)/bsxfun(@power,10,parameters_standard.MS.par(n_thetaDyn+(1),1))));
    maxErk_s(i) = max(max(bestFit_s(:,2,i)/bsxfun(@power,10,parameters_standard.MS.par(n_thetaDyn+(2),1)))); 
    
    %hierarchical approach
    sol = simulate_RafMekErk_relative(tsim,parameters_hierarchical.MS.par(:,1),kappa(i,:),[],options_simu);
    
    bestFit_h(:,:,i) = sol.y;
    
    maxMek_h(i) = max(max(bestFit_h(:,1,i)));
    maxErk_h(i) = max(max(bestFit_h(:,2,i)));
end

maxMekErk_s = repmat([max(maxMek_s) max(maxErk_s)],1,2);  
maxMekErk_h = repmat([max(maxMek_h) max(maxErk_h)],1,2); 

D_s.D.measurement{1}(:,[3,4,7,8]) = NaN;

figure(1)
for i = 1:n_kappa
    subplot(2,3,i)  
    
    if i == 1
        plot(D_h.D(i).t,D_h.D(i).my(:,1,1)/(c(:,1,1)*maxMekErk_h(1)),'og','Markersize',3);hold on
        plot(D_h.D(i).t,D_h.D(i).my(:,1,3)/(c(:,1,3)*maxMekErk_h(1)),'+g','Markersize',3);hold on
    else 
        plot(D_h.D(i).t,D_h.D(i).my(:,1,2)/(c(:,1,2)*maxMekErk_h(1)),'>g','Markersize',3);hold on
        plot(D_h.D(i).t,D_h.D(i).my(:,1,4)/(c(:,1,4)*maxMekErk_h(1)),'xg','Markersize',3);hold on
    end
    
    plot(tsim,bestFit_h(:,1,i)/maxMekErk_h(1),'Color',col_h,'LineWidth',2);hold on
     
    if i == 1
        plot(D_s.D.t{i},D_s.D.measurement{i}(:,1)/(bsxfun(@power,10,parameters_standard.MS.par(13,1))*maxMekErk_s(1)),'ok','Markersize',3);hold on
        plot(D_s.D.t{i},D_s.D.measurement{i}(:,5)/(bsxfun(@power,10,parameters_standard.MS.par(17,1))*maxMekErk_s(1)),'+k','Markersize',3);hold on
    else
        plot(D_s.D.t{i},D_s.D.measurement{i}(:,3)/(bsxfun(@power,10,parameters_standard.MS.par(15,1))*maxMekErk_s(1)),'>k','Markersize',3);hold on
        plot(D_s.D.t{i},D_s.D.measurement{i}(:,7)/(bsxfun(@power,10,parameters_standard.MS.par(19,1))*maxMekErk_s(1)),'xk','Markersize',3);hold on
    end
    
    plot(tsim,bestFit_s(:,1,i)/(bsxfun(@power,10,parameters_standard.MS.par(13,1))*maxMekErk_s(1)),'--','Color',col_s,'LineWidth',2);hold on
  
    title(title_vec{i},'FontSize',8,'FontWeight',setLay.lab.weight);
    xlabel(xlab{i},'FontSize',setLay.ax.size)
    ylabel(ylab{i},'FontSize',setLay.ax.size)
%     ylim([-0.2 1.5])
%     xlim([-0.1 10.1])
    
    axes = gca;
    if(i==1)
        set(axes,'Xtick',[0,2,4,6,8,10],'Ytick',[0,0.5,1,1.5],'XtickLabel',{''},'FontSize',14); 
    else        
        set(axes,'Xtick',[0,2,4,6,8,10],'Ytick',[0,0.5,1,1.5],'XtickLabel',{''},'YtickLabel',{''},'FontSize',14); 
    end
    
    subplot(2,3,i+3)
    
    if i == 1
        plot(D_h.D(i).t,D_h.D(i).my(:,2,1)/(c(:,2,1)*maxMekErk_h(2)),'og','Markersize',3);hold on
        plot(D_h.D(i).t,D_h.D(i).my(:,2,3)/(c(:,2,3)*maxMekErk_h(2)),'+g','Markersize',3);hold on
    else
        plot(D_h.D(i).t,D_h.D(i).my(:,2,2)/(c(:,2,2)*maxMekErk_h(2)),'>g','Markersize',3);hold on
        plot(D_h.D(i).t,D_h.D(i).my(:,2,4)/(c(:,2,4)*maxMekErk_h(2)),'xg','Markersize',3);hold on
    end
    
    plot(tsim,bestFit_h(:,2,i)/maxMekErk_h(2),'Color',col_h,'LineWidth',2);hold on
    
    if i == 1
        plot(D_s.D.t{i},D_s.D.measurement{i}(:,2)/(bsxfun(@power,10,parameters_standard.MS.par(14,1))*maxMekErk_s(2)),'ok','Markersize',3);hold on
        plot(D_s.D.t{i},D_s.D.measurement{i}(:,6)/(bsxfun(@power,10,parameters_standard.MS.par(18,1))*maxMekErk_s(2)),'+k','Markersize',3);hold on
    else
        plot(D_s.D.t{i},D_s.D.measurement{i}(:,4)/(bsxfun(@power,10,parameters_standard.MS.par(16,1))*maxMekErk_s(2)),'>k','Markersize',3);hold on
        plot(D_s.D.t{i},D_s.D.measurement{i}(:,8)/(bsxfun(@power,10,parameters_standard.MS.par(20,1))*maxMekErk_s(2)),'xk','Markersize',3);hold on
    end
    
    plot(tsim,bestFit_s(:,2,i)/(bsxfun(@power,10,parameters_standard.MS.par(14,1))*maxMekErk_s(2)),'--','Color',col_s,'LineWidth',2);
    
    xlabel(xlab{i+3},'FontSize',setLay.ax.size)
    ylabel(ylab{i+3},'FontSize',setLay.ax.size)
    ylim([-0.2 1.5])
    xlim([-0.1 10.1])
    
    axes = gca;
    if(i==1)
        set(axes,'Xtick',[0,2,4,6,8,10],'Ytick',[0,0.5,1,1.5],'FontSize',14); 
    else        
        set(axes,'Xtick',[0,2,4,6,8,10],'Ytick',[0,0.5,1,1.5],'YtickLabel',{''},'FontSize',14); 
    end
end

% % x0=1;
% % y0=1;
% % width=500;
% % height=600;
% % set(gcf,'units','points','position',[x0,y0,width,height]);
% 
%Convergence
starts_h = parameters_hierarchical.MS.logPost;
starts_h(abs(starts_h-starts_h(1))>1e-3) = NaN;
starts_h(~isnan(starts_h)) = 1;
n_conv_starts_s = nansum(starts_h);
c = figure(2);
plot(-parameters_hierarchical.MS.logPost,'o-','Color',col_h);
hold on
plot(-parameters_standard.MS.logPost,'o-','Color',col_s);
hold on
plot((1:2:n_conv_starts_s),-parameters_hierarchical.MS.logPost(1:2:n_conv_starts_s),'o-','Color',col_h);
xlim([0,300]);
%ylim([45,100]);
xlabel('optimizer runs','fontsize',setLay.ax.size)%,'FontName',setLay.lab.font)
ylabel('negative log-likelihood','fontsize',setLay.ax.size)%,'FontName',setLay.lab.font)
legend('hierarchical approach','standard approach');
set(gca,'FontSize',setLay.ax.size)%,'FontName',setLay.lab.font)
%title('Convergence','fontsize',setLay.tit.size)%,'FontName',setLay.lab.font);
pos = get(c,'Position');
%set(c,'Position',[pos(1) pos(2) 0.6*pos(3) 0.6*pos(4)]); 
set(c,'Position',[pos(1) pos(2) pos(3) pos(4)]);
% %print('-depsc2','-r1200','Convergence');
% x0=1;
% y0=1;
% width=500;
% height=600;
% set(gcf,'units','points','position',[x0,y0,width,height]);

% %Histogram cpu time
% cpu_time_n = parameters_standard.MS.t_cpu;
% cpu_time_n(cpu_time_n==0) = NaN;
% cpu_time_a = parameters_hierarchical.MS.t_cpu;
% cpu_time_a(cpu_time_a==0) = NaN;
% b = figure(3);
% %subplot(2,2,3)
% int = linspace(-1,5,100);
% intc = 0.5*(int(1:end-1)+int(2:end));
% h1 = histc(log10(cpu_time_n+1e-1),int);
% h1 = h1(1:end-1);
% h2 = histc(log10(cpu_time_a+1e-1),int);
% h2 = h2(1:end-1);
% bar(intc,h1,'BarWidth',1.0,'FaceColor',col_s);
% hold on
% bar(intc,h2,'BarWidth',0.7,'FaceColor',col_h); 
% %xlim(int([1,end]));
% xlabel('log_{10}(cpu time single start)[s]','fontsize',setLay.ax.size)%,'FontName',setLay.lab.font);
% ylabel('frequency','fontsize',setLay.ax.size)%,'FontName',setLay.lab.font);
% legend('standard approach','hierarchical approach');
% title('CPU time per optimization','fontsize',setLay.tit.size)%,'FontName',setLay.lab.font);
% set(gca,'FontSize',setLay.ax.size)%,'FontName',setLay.lab.font);
% %set(gca,'XTickLabel',{'10^{-1}' '10^0' '10^1' '10^2' '10^3' '10^4' '10^5'});
% pos = get(b,'Position');
% set(b,'Position',[pos(1) pos(2) 0.6*pos(3) pos(4)]); 
% %print('-depsc2','-r1200','Histogram');

% Boxplots cpu time
approach = {'standard','hierarchical'};
cpu_time_n = parameters_standard.MS.t_cpu;
cpu_time_n(cpu_time_n==0) = NaN;
cpu_time_a = parameters_hierarchical.MS.t_cpu;
cpu_time_a(cpu_time_a==0) = NaN;
bp = figure(4);
%boxplot([log10(cpu_time_n),log10(cpu_time_a)],approach,'Color', colors_sh,'OutlierSize',3,'Symbol','o');
boxplot([log10(cpu_time_n),log10(cpu_time_a)],approach,'OutlierSize',8,'Symbol','k.');hold on;
boxes = findobj(gca,'Tag','Box');
median = findobj(gca,'Tag','Median');
patch(get(boxes(1),'XData'),get(boxes(1),'YData'),col_h,'FaceAlpha',1); hold on;
patch(get(boxes(2),'XData'),get(boxes(2),'YData'),col_s,'FaceAlpha',1); hold on;
plot(get(median(1),'XData'),get(median(1),'YData'),'k-','Linewidth',0.5); hold on;
plot(get(median(2),'XData'),get(median(2),'YData'),'k-','Linewidth',0.5);
%xlabel('approach','fontsize',setLay.ax.size,'FontName',setLay.lab.font)
ylabel('cpu time single start [s]','fontsize',setLay.ax.size)%,'FontName',setLay.lab.font);
%title('CPU time per optimization','fontsize',setLay.tit.size)%,'FontName',setLay.lab.font);
set(gca,'YTick',(-1:1:3));
set(gca,'XTickLabel',{'standard','hierarchical'},'YTickLabel',{'10^{-1}' '10^0' '10^1' '10^2' '10^3'},'TickLabelInterpreter', 'tex');
set(gca,'FontSize',setLay.ax.size)%,'FontName',setLay.lab.font);
%set(gca,'YTickLabel',{'10^{-1}' '10^0' '10^1' '10^2' '10^3' '10^4' '10^5'},'TickLabelInterpreter', 'tex');
pos = get(bp,'Position');
%set(bp,'Position',[pos(1) pos(2) 0.6*pos(3) pos(4)]); 
set(bp,'Position',[pos(1) pos(2) pos(3) pos(4)]); 
% print('-depsc2','-r1200','Boxplot');
% x0=1;
% y0=1;
% width=500;
% height=600;
% set(gcf,'units','points','position',[x0,y0,width,height]);

% % % %Violinplot
% % % vp = figure(5);
% % % distributionPlot([log10(cpu_time_n),log10(cpu_time_a)],'color',color_vp,'histOpt',1); % histOpt=2 works better for uniform distributions than the default
% % % %xlabel('approach','fontsize',setLay.ax.size,'FontName',setLay.lab.font)
% % % %ylabel('log_{10}(cpu time single start)[s]','fontsize',setLay.ax.size)%,'FontName',setLay.lab.font);
% % % ylabel('cpu time single start [s]','fontsize',setLay.ax.size)%,'FontName',setLay.lab.font);
% % % title('CPU time per optimization','fontsize',setLay.tit.size)%,'FontName',setLay.lab.font);
% % % set(gca,'XTickLabel',{'standard' 'hierarchical'},'FontSize',setLay.ax.size)%,'FontName',setLay.lab.font);
% % % set(gca,'YTickLabel',{'10^{-1}' '10^0' '10^1' '10^2' '10^3' '10^4' '10^5'});%,'TickLabelInterpreter', 'tex');
% % % pos = get(vp,'Position');
% % % set(vp,'Position',[pos(1) pos(2) 0.6*pos(3) pos(4)]); 
% % % %print('-depsc2','-r1200','Violinplot');

%Cpu time/converged start
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
v = figure(6);
%subplot(2,2,4)
b_n = bar(vs_n,'BarWidth',0.5,'FaceColor',col_s);
hold on
b_a = bar(vs_a,'BarWidth',0.5,'FaceColor',col_h);
%set(gca, 'XTickLabel',{'standard','hierarchical'},'FontSize',setLay.ax.size)%,'FontName',setLay.lab.font);
set(gca,'YTick',(0:1:3));
set(gca,'FontSize',setLay.ax.size);
set(gca, 'XTickLabel',{'standard','hierarchical'},'YTickLabel',{'10^0' '10^1' '10^2' '10^3'},'TickLabelInterpreter', 'tex');
% legend(['total cpu time standard: ' num2str(cpu_time_n) 's'],...
%        ['total cpu time hierarchical: ' num2str(cpu_time_a) 's']);
%ylabel({'cpu time [s] per'; 'converged starts'},'fontsize',setLay.ax.size,'FontName',setLay.lab.font);
ylabel('cpu time [s] per converged start','fontsize',setLay.ax.size)%,'FontName',setLay.lab.font);
%title('CPU time per converged start','fontsize',setLay.tit.size)%,'FontName',setLay.lab.font);
pos = get(v,'Position');
%set(v,'Position',[pos(1) pos(2) 0.6*pos(3) 0.6*pos(4)]); 
set(v,'Position',[pos(1) pos(2) pos(3) pos(4)]); 
%print('-depsc2','-r1200','Retreat_CpuTime_Comparison');
% x0=1;
% y0=1;
% width=500;
% height=600;
% set(gcf,'units','points','position',[x0,y0,width,height]);
