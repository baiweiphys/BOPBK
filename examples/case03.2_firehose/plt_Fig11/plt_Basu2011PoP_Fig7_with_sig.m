% plt_Basu2011PoP_Fig7_wi.m

close all;
clear;
% clc;

%%
load('./PBK_Fig7_data/Basu2011PoPFig7_TizTix=5_km3.mat');
k_IC_TizTix5_km3 = pas;
w_IC_TizTix5_km3 = wws;
% wr_IC_TizTix5_km3 = real(w_IC_TizTix5_km3);
wi_IC_TizTix5_km3 = imag(w_IC_TizTix5_km3);
%
load('./PBK_Fig7_data/Basu2011PoPFig7_TizTix=5_sig=0.5_km3.mat');
k_IC_TizTix5_sig0dot5_km3 = pas;
w_IC_TizTix5_sig0dot5_km3 = wws;
% wr_IC_TizTix5_sig0dot5_km3 = real(w_IC_TizTix5_sig0dot5_km3);
wi_IC_TizTix5_sig0dot5_km3 = imag(w_IC_TizTix5_sig0dot5_km3);
%
load('./PBK_Fig7_data/Basu2011PoPFig7_TizTix=5_sig=1_km3.mat');
k_IC_TizTix5_sig1_km3 = pas;
w_IC_TizTix5_sig1_km3 = wws;
% wr_IC_TizTix5_sig1_km3 = real(w_IC_TizTix5_sig1_km3);
wi_IC_TizTix5_sig1_km3 = imag(w_IC_TizTix5_sig1_km3);

%%
load('./PBK_Fig7_data/Basu2011PoPFig7_TizTix=5_bm.mat');
k_IC_TizTix5_bm = pas;
w_IC_TizTix5_bm = wws;
% wr_IC_TizTix5_bm = real(w_IC_TizTix5_bm);
wi_IC_TizTix5_bm = imag(w_IC_TizTix5_bm);
%
load('./PBK_Fig7_data/Basu2011PoPFig7_TizTix=5_sig=0.5_bm.mat');
k_IC_TizTix5_sig0dot5_bm = pas;
w_IC_TizTix5_sig0dot5_bm = wws;
% wr_IC_TizTix5_sig0dot5_bm = real(w_IC_TizTix5_sig0dot5_bm);
wi_IC_TizTix5_sig0dot5_bm = imag(w_IC_TizTix5_sig0dot5_bm);
%
load('./PBK_Fig7_data/Basu2011PoPFig7_TizTix=5_sig=1_bm.mat');
k_IC_TizTix5_sig1_bm = pas;
w_IC_TizTix5_sig1_bm = wws;
% wr_IC_TizTix5_sig1_bm = real(w_IC_TizTix5_sig1_bm);
wi_IC_TizTix5_sig1_bm = imag(w_IC_TizTix5_sig1_bm);

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
h=figure('unit','normalized','Position',[0.01 0.45 0.43 0.45],'DefaultAxesFontSize',24);
ha = tight_subplot(1,1,[0.03 .03],[0.17 0.03],[0.15 0.04]);
%
axes(ha(1)); 
plot(k_IC_TizTix5_km3,wi_IC_TizTix5_km3(:,1,1),'-','Color',pltc(8,:),'linewidth',3);
hold on
plot(k_IC_TizTix5_sig0dot5_km3,wi_IC_TizTix5_sig0dot5_km3(:,1,1),'-','Color',pltc(2,:),'linewidth',3);
plot(k_IC_TizTix5_sig1_km3,wi_IC_TizTix5_sig1_km3(:,1,1),'-','Color',pltc(3,:),'linewidth',3);
%
plot(k_IC_TizTix5_bm,wi_IC_TizTix5_bm(:,1,1),'-.','Color',pltc(8,:),'linewidth',3);
plot(k_IC_TizTix5_sig0dot5_bm,wi_IC_TizTix5_sig0dot5_bm(:,1,1),'-.','Color',pltc(2,:),'linewidth',3);
plot(k_IC_TizTix5_sig1_bm,wi_IC_TizTix5_sig1_bm(:,1,1),'-.','Color',pltc(3,:),'linewidth',3);
%
grid on;
xlabel('$k_{\perp} \rho_{i}$','Interpreter','latex','FontSize',32);
ylabel('$\omega_i/\omega_{ci}$','Interpreter','latex','FontSize',32);
%
legend('PBK($\kappa_{\parallel(e,i)}=3$, $\kappa_{\perp(e,i)}=200$, $\sigma_i=0$)',...
       'PBK($\kappa_{\parallel(e,i)}=3$, $\kappa_{\perp(e,i)}=200$, $\sigma_i=0.5$)',...
       'PBK($\kappa_{\parallel(e,i)}=3$, $\kappa_{\perp(e,i)}=200$, $\sigma_i=1$)',...
       'BM($\kappa_{\parallel(e,i)}=\infty$, $\kappa_{\perp(e,i)}=\infty$, $\sigma_i=0$)',...
       'BM($\kappa_{\parallel(e,i)}=\infty$, $\kappa_{\perp(e,i)}=\infty$, $\sigma_i=0.5$)',...
       'BM($\kappa_{\parallel(e,i)}=\infty$, $\kappa_{\perp(e,i)}=\infty$, $\sigma_i=1$)',...
       'Interpreter','latex','FontSize',20,'Location','best');
legend('boxoff');
ylim([0,0.3]);
xlim([0 2]);
% ylim([-5e-3,1e-1]);

set_XY_Tick;

%% save figure
savefig('benchmark_Basu2011PoP_fig7_sig.fig');

set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],...
  'PaperSize',[screenposition(3:4)]);

% print(gcf,['benchmark_Basu2011PoP_fig7','.eps'],'-depsc');
print(gcf,['benchmark_Basu2011PoP_fig7_sig','.pdf'],'-dpdf');
