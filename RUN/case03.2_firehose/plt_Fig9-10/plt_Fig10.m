% plt_Basu2011PoP_Fig7_wi.m

close all;
clear;
% clc;

%%
load('./PBK_Fig7_data/Basu2011PoPFig7_TizTix=5_km3.mat');
k_IC_TizTix5_km3 = pas;
w_IC_TizTix5_km3 = wws;
wr_IC_TizTix5_km3 = real(w_IC_TizTix5_km3);
wi_IC_TizTix5_km3 = imag(w_IC_TizTix5_km3);

%
load('./PBK_Fig7_data/Basu2011PoPFig7_TizTix=5_bm.mat');
k_IC_TizTix5_bm = pas;
w_IC_TizTix5_bm = wws;
wr_IC_TizTix5_bm = real(w_IC_TizTix5_bm);
wi_IC_TizTix5_bm = imag(w_IC_TizTix5_bm);
%
%
load('./PBK_Fig7_data/Basu2011PoPFig7_TizTix=2_km3.mat');
k_IC_TizTix2_km3 = pas;
w_IC_TizTix2_km3 = wws;
wr_IC_TizTix2_km3 = real(w_IC_TizTix2_km3);
wi_IC_TizTix2_km3 = imag(w_IC_TizTix2_km3);
%
load('./PBK_Fig7_data/Basu2011PoPFig7_TizTix=2_bm.mat');
k_IC_TizTix2_bm = pas;
w_IC_TizTix2_bm = wws;
wr_IC_TizTix2_bm = real(w_IC_TizTix2_bm);
wi_IC_TizTix2_bm = imag(w_IC_TizTix2_bm);



%%
% load csv file
%
Basu2011PoP_Fig7_wi_TizTix5 = csvread('./Basu2011PoP_Fig7_data/Basu2011PoP_Fig7_wi_TizTix5.csv',1,0);
k_wi_Basu2011PoP_Fig7_TizTix5 = Basu2011PoP_Fig7_wi_TizTix5(:,1);
wi_Basu2011PoP_Fig7_TizTix5_km3 = Basu2011PoP_Fig7_wi_TizTix5(:,2);
wi_Basu2011PoP_Fig7_TizTix5_bm = Basu2011PoP_Fig7_wi_TizTix5(:,3);
%
Basu2011PoP_Fig7_wi_TizTix2 = csvread('./Basu2011PoP_Fig7_data/Basu2011PoP_Fig7_wi_TizTix2.csv',1,0);
k_wi_Basu2011PoP_Fig7_TizTix2 = Basu2011PoP_Fig7_wi_TizTix2(:,1);
wi_Basu2011PoP_Fig7_TizTix2_km3 = Basu2011PoP_Fig7_wi_TizTix2(:,2);
wi_Basu2011PoP_Fig7_TizTix2_bm = Basu2011PoP_Fig7_wi_TizTix2(:,3);
%
%
Basu2011PoP_Fig7_wr_TizTix5 = csvread('./Basu2011PoP_Fig7_data/Basu2011PoP_Fig7_wr_TizTix5.csv',1,0);
k_wr_Basu2011PoP_Fig7_TizTix5 = Basu2011PoP_Fig7_wr_TizTix5(:,1);
wr_Basu2011PoP_Fig7_TizTix5_km3 = Basu2011PoP_Fig7_wr_TizTix5(:,2);
wr_Basu2011PoP_Fig7_TizTix5_bm = Basu2011PoP_Fig7_wr_TizTix5(:,3);

%%
Basu2011PoP_Fig7_wr_TizTix2 = csvread('./Basu2011PoP_Fig7_data/Basu2011PoP_Fig7_wr_TizTix2.csv',1,0);
k_wr_Basu2011PoP_Fig7_TizTix2 = Basu2011PoP_Fig7_wr_TizTix2(:,1);
wr_Basu2011PoP_Fig7_TizTix2_km3 = Basu2011PoP_Fig7_wr_TizTix2(:,2);
wr_Basu2011PoP_Fig7_TizTix2_bm = Basu2011PoP_Fig7_wr_TizTix2(:,3);

%%
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
h=figure('unit','normalized','Position',[0.01 0.45 0.37 0.6],'DefaultAxesFontSize',24);
ha = tight_subplot(2,1,[0.03 .03],[0.12 0.03],[0.19 0.05]);
% subplot(2,1,1);
axes(ha(1)); 
plot(k_IC_TizTix2_km3,wr_IC_TizTix2_km3(:,1,1),'-','Color',pltc(1,:),'linewidth',3);
hold on
plot(k_IC_TizTix2_bm,wr_IC_TizTix2_bm(:,1,1),'-','Color',pltc(4,:),'linewidth',3);
plot(k_IC_TizTix5_km3,wr_IC_TizTix5_km3(:,1,1),'-','Color',pltc(3,:),'linewidth',3);
plot(k_IC_TizTix5_bm,wr_IC_TizTix5_bm(:,1,1),'-','Color',pltc(2,:),'linewidth',3);

%
plot(k_wr_Basu2011PoP_Fig7_TizTix5,wr_Basu2011PoP_Fig7_TizTix5_km3,'--','Color',pltc(8,:),'LineWidth',3);
plot(k_wr_Basu2011PoP_Fig7_TizTix5,wr_Basu2011PoP_Fig7_TizTix5_bm,'--','Color',pltc(8,:),'LineWidth',3);
plot(k_wr_Basu2011PoP_Fig7_TizTix2,wr_Basu2011PoP_Fig7_TizTix2_km3,'--','Color',pltc(8,:),'LineWidth',3);
plot(k_wr_Basu2011PoP_Fig7_TizTix2,wr_Basu2011PoP_Fig7_TizTix2_bm,'--','Color',pltc(8,:),'LineWidth',3);
%
% xlabel('ck/\omega_{pe}');
ylabel('$\omega_r/\omega_{ci}$','Interpreter','latex');
% title('OFHI: $\theta=45^\circ$','Interpreter','latex','FontSize',30);
%
xlim([0 3]);
ylim([1.15,1.6]);

set(ha(1),'XTickLabel',''); 
set_XY_Tick;

%%
% subplot(2,1,2);
axes(ha(2)); 
plot(k_IC_TizTix2_km3,wi_IC_TizTix2_km3(:,1,1),'-','Color',pltc(1,:),'linewidth',3);
hold  on
plot(k_IC_TizTix2_bm,wi_IC_TizTix2_bm(:,1,1),'-','Color',pltc(4,:),'linewidth',3);
plot(k_IC_TizTix5_km3,wi_IC_TizTix5_km3(:,1,1),'-','Color',pltc(3,:),'linewidth',3);
plot(k_IC_TizTix5_bm,wi_IC_TizTix5_bm(:,1,1),'-','Color',pltc(2,:),'linewidth',3);

plot(k_wi_Basu2011PoP_Fig7_TizTix5,wi_Basu2011PoP_Fig7_TizTix5_km3,'--','Color',pltc(8,:),'LineWidth',3);
plot(k_wi_Basu2011PoP_Fig7_TizTix5,wi_Basu2011PoP_Fig7_TizTix5_bm,'--','Color',pltc(8,:),'LineWidth',3);
plot(k_wi_Basu2011PoP_Fig7_TizTix2,wi_Basu2011PoP_Fig7_TizTix2_km3,'--','Color',pltc(8,:),'LineWidth',3);
plot(k_wi_Basu2011PoP_Fig7_TizTix2,wi_Basu2011PoP_Fig7_TizTix2_bm,'--','Color',pltc(8,:),'LineWidth',3);

%
grid on;
xlabel('$k_{\perp} \rho_{i}$','Interpreter','latex');
ylabel('$\omega_i/\omega_{ci}$','Interpreter','latex');
%
legend('$T_{\parallel i}/T_{\perp i}=2$, PBK($\kappa_{\parallel(e,i)}=3$, $\kappa_{\perp(e,i)}=200$)',...
       '$T_{\parallel i}/T_{\perp i}=2$, BM($\kappa_{\parallel(e,i)}=\infty$, $\kappa_{\perp(e,i)}=\infty$)',...
       '$T_{\parallel i}/T_{\perp i}=5$, PBK($\kappa_{\parallel(e,i)}=3$, $\kappa_{\perp(e,i)}=200$)',...
       '$T_{\parallel i}/T_{\perp i}=5$, BM($\kappa_{\parallel(e,i)}=\infty$, $\kappa_{\perp(e,i)}=\infty$)',...
       'Basu(2011)',...
       'Interpreter','latex','FontSize',16.5,'Location','best');
legend('boxoff');
ylim([0,0.3]);
xlim([0 3]);
% ylim([-5e-3,1e-1]);

set_XY_Tick;

%% save figure
savefig('benchmark_Basu2011PoP_fig7.fig');

set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],...
  'PaperSize',[screenposition(3:4)]);

% print(gcf,['benchmark_Basu2011PoP_fig7','.eps'],'-depsc');
print(gcf,['benchmark_Basu2011PoP_fig7','.pdf'],'-dpdf');
