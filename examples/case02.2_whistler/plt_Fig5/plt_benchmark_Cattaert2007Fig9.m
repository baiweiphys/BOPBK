% plt_benchmark_Cattaert2007Fig9.m

clear all;
% clc;
close all;


load('./pbkData/saveData_Cattaert2007Fig9_TexTez0dot5.mat');
kk_km2_TexTez0dot5 = pas;
wr_km2_TexTez0dot5 = real(wws);
wi_km2_TexTez0dot5 = imag(wws);

load('./pbkData/saveData_Cattaert2007Fig9_TexTez1.mat');
kk_km2_TexTez1 = pas;
wr_km2_TexTez1 = real(wws);
wi_km2_TexTez1 = imag(wws);

load('./pbkData/saveData_Cattaert2007Fig9_TxTz2.mat');
kk_km2_TxTz2 = pas;
wr_km2_TxTz2 = real(wws);
wi_km2_TxTz2 = imag(wws);

load('./pbkData/saveData_Cattaert2007Fig9_TxTz3.mat');
kk_km2_TxTz3 = pas;
wr_km2_TxTz3 = real(wws);
wi_km2_TxTz3 = imag(wws);


% load csv file
wi = csvread('./Cattaert2007Fig9_data/wi_Cattaert2007Fig9.csv',1,0);
wr = csvread('./Cattaert2007Fig9_data/wr_Cattaert2007Fig9.csv',1,0);

wr_x = wr(:,1);
wr_TxTz0dot5 = wr(:,2);
wr_TxTz1 = wr(:,3);
wr_TxTz2 = wr(:,4);
wr_TxTz3 = wr(:,5);

wi_x = wi(:,1);
wi_TxTz0dot5 = wi(:,2);
wi_TxTz1 = wi(:,3);
wi_TxTz2 = wi(:,4);
wi_TxTz3 = wi(:,5);


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
rootNo_TxTz0dot5 = 36:37;
rootNo_TxTz1  = 36:37;
rootNo_TxTz2  = 34:36;
rootNo_TxTz3  = 2:30;

elem = 2;

h=figure('unit','normalized','Position',[0.01 0.45 0.35 0.57],'DefaultAxesFontSize',22);
% subplot(2,1,1);
ha = tight_subplot(2,1,[0.03 .03],[0.13 .02],[0.19 0.05]);
axes(ha(1)); 
plot(kk_km2_TexTez0dot5,abs(wr_km2_TexTez0dot5),'-','Color',pltc(1,:),'linewidth',3);
hold on;
plot(kk_km2_TexTez1,abs(wr_km2_TexTez1),'-','Color',pltc(2,:),'linewidth',3);
plot(kk_km2_TxTz2,abs(wr_km2_TxTz2),'-','Color',pltc(3,:),'linewidth',3);
plot(kk_km2_TxTz3,abs(wr_km2_TxTz3),'-','Color',pltc(4,:),'linewidth',3);
%
plot(wr_x,wr_TxTz0dot5,'k--','LineWidth',2);
plot(wr_x,wr_TxTz1,'k--','LineWidth',2);
plot(wr_x,wr_TxTz2,'k--','LineWidth',2);
plot(wr_x,wr_TxTz3,'k--','LineWidth',2);
ylabel('$\omega_r/|\omega_{ce}|$','Interpreter','latex','FontSize',32);
% xlabel('k\rho_{ce}');
xlim([0 1]);
ylim([1e-3 0.9]);
legend('$T_{\perp e}/T_{\parallel e}=0.5$','$T_{\perp e}/T_{\parallel e}=1.0$',...
    '$T_{\perp e}/T_{\parallel e}=2.0$','$T_{\perp e}/T_{\parallel e}=3.0$',...
    'Cattaert(2007)','Interpreter','latex','FontSize',22,'Location','southeast');
legend('boxoff');
set(ha(1),'XTickLabel',''); 
set_XY_Tick;

% subplot(2,1,2);
axes(ha(2)); 
plot(kk_km2_TexTez0dot5,wi_km2_TexTez0dot5,'-','Color',pltc(1,:),'linewidth',3);
hold on;
plot(kk_km2_TexTez1,wi_km2_TexTez1,'-','Color',pltc(2,:),'linewidth',3);
plot(kk_km2_TxTz2,wi_km2_TxTz2,'-','Color',pltc(3,:),'linewidth',3);
plot(kk_km2_TxTz3,wi_km2_TxTz3,'-','Color',pltc(4,:),'linewidth',3);
%
plot(wi_x,wi_TxTz0dot5,'k--','LineWidth',2);
plot(wi_x,wi_TxTz1,'k--','LineWidth',2);
plot(wi_x(1:end-1),wi_TxTz2(1:end-1),'k--','LineWidth',2);
plot(wi_x,wi_TxTz3,'k--','LineWidth',2);
pl = line([0,1],[0,0]);
pl.Color = 'black';
pl.LineStyle = '--';
xlabel('$k\rho_{ce}$','Interpreter','latex','FontSize',32);
ylabel('$\omega_i/|\omega_{ce}|$','Interpreter','latex','FontSize',32);
xlim([0 1]);
ylim([-0.04,0.03]);
set_XY_Tick;

%%
% save figure
savefig('benchmark_Cattaert2007Fig9.fig');

set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],...
  'PaperSize',[screenposition(3:4)]);

% print(gcf,['benchmark_Cattaert2007Fig9','.eps'],'-depsc');
print(gcf,['benchmark_Cattaert2007Fig9','.pdf'],'-dpdf');


