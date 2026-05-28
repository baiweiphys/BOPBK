% plt_benchmark_Cattaert2007Fig9.m

clear all;
% clc;
close all;


% load PBKi20 
load('./pbkData/Santos2015Fig6a_IC_TixTiz=7_pbki20.mat');
kk_TixTiz7_pbk20 = pas;
wr_TixTiz7_pbk20 = real(wws);
wi_TixTiz7_pbk20 = imag(wws);
%
load('./pbkData/Santos2015Fig6a_IC_TixTiz=5_pbki20.mat');
kk_TixTiz5_pbk20 = pas;
wr_TixTiz5_pbk20 = real(wws);
wi_TixTiz5_pbk20 = imag(wws);
%
load('./pbkData/Santos2015Fig6a_IC_TixTiz=3_pbki20.mat');
kk_TixTiz3_pbk20 = pas;
wr_TixTiz3_pbk20 = real(wws);
wi_TixTiz3_pbk20 = imag(wws);

% load PBKi5 
load('./pbkData/Santos2015Fig6b_IC_TixTiz=7_pbki5.mat');
kk_TixTiz7_pbk5 = pas;
wr_TixTiz7_pbk5 = real(wws);
wi_TixTiz7_pbk5 = imag(wws);
%
load('./pbkData/Santos2015Fig6b_IC_TixTiz=5_pbki5.mat');
kk_TixTiz5_pbk5 = pas;
wr_TixTiz5_pbk5 = real(wws);
wi_TixTiz5_pbk5 = imag(wws);
%
load('./pbkData/Santos2015Fig6b_IC_TixTiz=3_pbki5.mat');
kk_TixTiz3_pbk5 = pas;
wr_TixTiz3_pbk5 = real(wws);
wi_TixTiz3_pbk5 = imag(wws);



% load csv file
wi_pbk20 = csvread('./Santos2015PoP_data/Santos2015PoP_Fig6a.csv',1,0);
x_wi_fig6a = wi_pbk20(:,1);
wi_TixTiz7_fig6a = wi_pbk20(:,2);
wi_TixTiz5_fig6a = wi_pbk20(:,3);
wi_TixTiz3_fig6a = wi_pbk20(:,4);

%
wi_pbk5 = csvread('./Santos2015PoP_data/Santos2015PoP_Fig6b.csv',1,0);
x_wi_fig6b = wi_pbk5(:,1);
wi_TixTiz7_fig6b = wi_pbk5(:,2);
wi_TixTiz5_fig6b = wi_pbk5(:,3);
wi_TixTiz3_fig6b = wi_pbk5(:,4);



% for plots with different colors
pltc=[0.0  1.0  0.0
  1.0  0.0  0.0
  0.2  0.2  1.0
  0.8  0.8  0.0
  1.0  0.6  0.0
  0.9  0.0  0.9
  0.0  0.8  0.8
  0.0  0.0  0.0
  0.6  0.0  0.0
  0.4  0.7  0.4 
  0.0  0.0  0.5 
  0.6  0.0  0.6 
  0.0  0.5  1.0
  ];

%% plots

% elem = 2;

% h=figure('unit','normalized','Position',[0.01 0.45 0.35 0.57],'DefaultAxesFontSize',22);
% ha = tight_subplot(2,1,[0.03 .03],[0.13 .02],[0.19 0.05]);

% h=figure('unit','normalized','Position',[0.02 0.45 0.55 0.45],'DefaultAxesFontSize',23);
% ha = tight_subplot(1,2,[0.0 0.02],[0.16 0.09],[0.08 0.1]);

h=figure('unit','normalized','Position',[0.02 0.45 0.55 0.46],'DefaultAxesFontSize',24);
ha = tight_subplot(1,2,[0.0 0.02],[0.16 0.09],[0.1 0.05]);
%% subplot(2,1,1);
axes(ha(1)); 
plot(kk_TixTiz7_pbk20,wi_TixTiz7_pbk20,'-','Color',pltc(1,:),'linewidth',4);
hold on;
plot(kk_TixTiz5_pbk20,wi_TixTiz5_pbk20,'-','Color',pltc(2,:),'linewidth',4);
plot(kk_TixTiz3_pbk20,wi_TixTiz3_pbk20,'-','Color',pltc(4,:),'linewidth',4);
%
plot(x_wi_fig6a,wi_TixTiz7_fig6a,'k--','LineWidth',4);
plot(x_wi_fig6a,wi_TixTiz5_fig6a,'k--','LineWidth',4);
plot(x_wi_fig6a,wi_TixTiz3_fig6a,'k--','LineWidth',4);
% pl = line([0,1],[0,0]);
% pl.Color = 'black';
% pl.LineStyle = '--';
xlabel('$k v_A/\omega_{ci}$','Interpreter','latex','FontSize',34);
ylabel('$\omega_i/\omega_{ci}$','Interpreter','latex','FontSize',34);
legend('$T_{\perp i}/T_{\parallel i}=7$','$T_{\perp i}/T_{\parallel i}=5$',...
    '$T_{\perp i}/T_{\parallel i}=3$',...
    'Santos(2015)','Interpreter','latex','FontSize',24,'Location','northeast');
legend('boxoff');
xlim([0 4]);
ylim([0.0,1.2]);
set_XY_Tick;
title(['(a) $\kappa_{\parallel i}=\kappa_{\perp i}=20$ ', '($\kappa_{\parallel e}=\kappa_{\perp e}=\infty$)'],...
       'Interpreter','latex','FontSize',25);

%% subplot(2,1,2);
axes(ha(2)); 
plot(kk_TixTiz7_pbk5,wi_TixTiz7_pbk5,'-','Color',pltc(1,:),'linewidth',4);
hold on;
plot(kk_TixTiz5_pbk5,wi_TixTiz5_pbk5,'-','Color',pltc(2,:),'linewidth',4);
plot(kk_TixTiz3_pbk5,wi_TixTiz3_pbk5,'-','Color',pltc(4,:),'linewidth',4);
%
plot(x_wi_fig6b,wi_TixTiz7_fig6b,'k--','LineWidth',4);
plot(x_wi_fig6b,wi_TixTiz5_fig6b,'k--','LineWidth',4);
plot(x_wi_fig6b,wi_TixTiz3_fig6b,'k--','LineWidth',4);
% pl = line([0,1],[0,0]);
% pl.Color = 'black';
% pl.LineStyle = '--';
xlabel('$k v_A/\omega_{ci}$','Interpreter','latex','FontSize',34);
% legend('$T_{\perp i}/T_{\parallel i}=7$','$T_{\perp i}/T_{\parallel i}=5$',...
%     '$T_{\perp i}/T_{\parallel i}=3$',...
%     'Santos(2015)','Interpreter','latex','FontSize',22,'Location','northeast');
% legend('boxoff');
xlim([0 4]);
ylim([0.0,1.2]);
set_XY_Tick;
set(ha(2),'YTickLabel',''); 
title(['(b) $\kappa_{\parallel i}=\kappa_{\perp i}=5$ ', '($\kappa_{\parallel e}=\kappa_{\perp e}=\infty$)'],...
       'Interpreter','latex','FontSize',25);

%%
% save figure
savefig('benchmark_Santos2015Fig6ab.fig');

set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],...
  'PaperSize',[screenposition(3:4)]);

% print(gcf,['benchmark_Santos2015Fig1','.eps'],'-depsc');
print(gcf,['benchmark_Santos2015Fig6ab','.pdf'],'-dpdf');


