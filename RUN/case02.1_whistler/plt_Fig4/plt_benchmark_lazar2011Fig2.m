% plot_benchmark_Cattaert2007Fig1.m

close all;
clear;
clc;

% load data
load('./BOPBK_data/saveData_Electron_Whistler_bm8_Lazar2011Fig2.mat');
w_bm = wws;
k0_bm = pas;
%
load('./BOPBK_data/saveData_Electron_Whistler_pbk2_Lazar2011Fig2.mat');
w_pbk2 = wws;
k0_pbk2 = pas;
%
load('./BOPBK_data/saveData_Electron_Whistler_pbk6_Lazar2011Fig2.mat');
w_pbk6 = wws;
k0_pbk6 = pas;



%% Load csv file
% load csv file of Fig2_bm
wr = csvread('./Lazar2011_data/Fig2_bm/wr_bm_Lazar2011Fig2.csv',1,0);
wi = csvread('./Lazar2011_data/Fig2_bm/wi_bm_Lazar2011Fig2.csv',1,0);
%
kr_bm_lazar2011Fig2 = wr(:,1);
wr_bm_lazar2011Fig2 = wr(:,2);
ki_bm_lazar2011Fig2 = wi(:,1);
wi_bm_lazar2011Fig2 = wi(:,2);

% load csv file of Fig2_pbk2
wr = csvread('./Lazar2011_data/Fig2_pbk3/wr_pbk3_Lazar2011Fig2.csv',1,0);
wi = csvread('./Lazar2011_data/Fig2_pbk3/wi_pbk3_Lazar2011Fig2.csv',1,0);
%
kr_pbk3_lazar2011Fig2 = wr(:,1);
wr_pbk3_lazar2011Fig2 = wr(:,2);
ki_pbk3_lazar2011Fig2 = wi(:,1);
wi_pbk3_lazar2011Fig2 = wi(:,2);

% load csv file of Fig2_pbk2
wr = csvread('./Lazar2011_data/Fig2_pbk7/wr_pbk7_Lazar2011Fig2.csv',1,0);
wi = csvread('./Lazar2011_data/Fig2_pbk7/wi_pbk7_Lazar2011Fig2.csv',1,0);
%
kr_pbk7_lazar2011Fig2 = wr(:,1);
wr_pbk7_lazar2011Fig2 = wr(:,2);
ki_pbk7_lazar2011Fig2 = wi(:,1);
wi_pbk7_lazar2011Fig2 = wi(:,2);


%
pltc=[0.0  1.0  0.0 % different colors
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
h=figure('unit','normalized','Position',[0.01 0.45 0.35 0.6],'DefaultAxesFontSize',20);
ha = tight_subplot(2,1,[0.03 .03],[0.12 0.03],[0.19 0.05]);
% subplot(2,1,1);
axes(ha(1)); 
% to benchmark the bm real part curve from lazar2011 Fi2 
plot(k0_bm,real(w_bm),'-','Color',pltc(1,:),'linewidth',3); 
hold on;
% to benchmark the pbk  real part curve from lazar2011 Fi2
plot(k0_pbk2,real(w_pbk2),'-','Color',pltc(2,:),'linewidth',3);  
plot(k0_pbk6,real(w_pbk6),'-','Color',pltc(3,:),'linewidth',3);  

%
plot(kr_bm_lazar2011Fig2,wr_bm_lazar2011Fig2,'k--','Markersize',12,'linewidth',3);
plot(kr_pbk3_lazar2011Fig2,wr_pbk3_lazar2011Fig2,'k--','Markersize',10,'linewidth',3);
plot(kr_pbk7_lazar2011Fig2(1:end-3),wr_pbk7_lazar2011Fig2(1:end-3),'k--','Markersize',10,'linewidth',3);
ylabel('$\omega_r/\omega_{pe}$','Fontsize',33,'Interpreter','latex');
% xlabel('k\rho_{ce}');
% xlim([1e-3 1e5]);
% ylim([0 3]);
legend('BO-PBK(BM, $J=8$)','BO-PBK(PBK, $\kappa_{\parallel,\perp(e,p)}=2$)', ...
       'BO-PBK(PBK, $\kappa_{\parallel,\perp(e,p)}=6$)', ...
       'Lazer(2010)',...
       'Interpreter','latex','FontSize',20,'Location','southeast');
legend('boxoff');
ylim([0.0,0.01]);
set(ha(1),'XTickLabel',''); 
set_XY_Tick;

% subplot(2,1,2);
axes(ha(2)); 
% to benchmark the bm image part curve from lazar2011 Fi2 
plot(k0_bm,imag(w_bm),'-','Color',pltc(1,:),'linewidth',3); % bm
hold on;
% to benchmark the pbk real part curve from lazar2011 Fi2 
plot(k0_pbk2,imag(w_pbk2),'-','Color',pltc(2,:),'linewidth',3); % pbk2
plot(k0_pbk6,imag(w_pbk6),'-','Color',pltc(3,:),'linewidth',3); % pbk7
%
plot(ki_bm_lazar2011Fig2,wi_bm_lazar2011Fig2,'k--','Markersize',10,'linewidth',3);
plot(ki_pbk3_lazar2011Fig2,wi_pbk3_lazar2011Fig2,'k--','Markersize',10,'linewidth',3);
plot(ki_pbk7_lazar2011Fig2(1:end-2),wi_pbk7_lazar2011Fig2(1:end-2),'k--','Markersize',10,'linewidth',3);
ylabel('$\omega_i/\omega_{pe}$','Fontsize',33,'Interpreter','latex');
xlabel('$ck/\omega_{pe}$','Fontsize',33,'Interpreter','latex');
% xlim([0,0.37]);
ylim([0.0,0.014]);
set_XY_Tick;


% save figure
savefig('benchmark_lazar2010Fig2.fig');

set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],...
  'PaperSize',[screenposition(3:4)]);

% print(gcf,['benchmark_lazar2010Fig2','.eps'],'-depsc');
print(gcf,['benchmark_lazar2010Fig2','.pdf'],'-dpdf');
