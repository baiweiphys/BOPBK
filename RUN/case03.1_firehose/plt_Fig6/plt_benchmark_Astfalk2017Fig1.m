%plt_benchmark_Astfalk2017Fig1

close all;
clear;
% clc;

%
load('./BOPBK_data/PFHI_bm/saveData_Astfalk2017Fig1_PFHI_bm.mat');
k_PFHI_bm = pas;
w_PFHI_bm = wws(:,:,1);
realw_PFHI_bm = real(w_PFHI_bm);
realw_PFHI_bm  = abs(realw_PFHI_bm);
imagw_PFHI_bm = imag(w_PFHI_bm);

%
load('./BOPBK_data/OFHI_th45_bm/saveData_Astfalk2017Fig1_OFHI_bm.mat');
k_OFHI_bm = pas;
w_OFHI_bm = wws;
realw_OFHI_bm = real(w_OFHI_bm);
imagw_OFHI_bm = imag(w_OFHI_bm);

%%
% load csv file
wr_PFHI_Astfalk2017Fig1 = csvread('./Astfalk2017Fig1_data/Astfalk2017Fig1_wr_th0_bm.csv',1,0);
wi_PFHI_Astfalk2017Fig1 = csvread('./Astfalk2017Fig1_data/Astfalk2017Fig1_wi_th0_bm.csv',1,0);
wi_OFHI_Astfalk2017Fig1 = csvread('./Astfalk2017Fig1_data/Astfalk2017Fig1_wi_th45_bm.csv',1,0);

k_wr_th0_Astfalk2017Fig1 = wr_PFHI_Astfalk2017Fig1(:,1);
wr_th0_Astfalk2017Fig1 = wr_PFHI_Astfalk2017Fig1(:,2);

k_wi_th0_Astfalk2017Fig1 = wi_PFHI_Astfalk2017Fig1(:,1);
wi_th0_Astfalk2017Fig1 = wi_PFHI_Astfalk2017Fig1(:,2);

k_wi_th45_Astfalk2017Fig1 = wi_OFHI_Astfalk2017Fig1(:,1);
wi_th45_Astfalk2017Fig1 = wi_OFHI_Astfalk2017Fig1(:,2);

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

h=figure('unit','normalized','Position',[0.01 0.45 0.35 0.54],'DefaultAxesFontSize',22);
ha = tight_subplot(2,1,[0.03 .03],[0.137 .02],[0.19 0.05]);
axes(ha(1)); 
plot(k_PFHI_bm,realw_PFHI_bm,'-','Color',pltc(1,:),'linewidth',3);
hold on;
plot(k_OFHI_bm,realw_OFHI_bm,'-','Color',pltc(2,:),'linewidth',3);
%
plot(k_wr_th0_Astfalk2017Fig1,wr_th0_Astfalk2017Fig1,'k--','LineWidth',3);
% xlabel('ck/\omega_{pe}');
ylabel('$\omega_r/\omega_{ci}$','Interpreter','latex','FontSize',32);

% xlim([0 3]);
ylim([0 1.2]);
legend('$\theta=0^\circ$, BM($J_{e,i}=8$)', ...
       '$\theta=45^\circ$, BM($J_{e,i}=8$)', ...
       'LEOPARD',...
       'Interpreter','latex','FontSize',22,'Location','northwest');
legend('boxoff');
set(ha(1),'XTickLabel',''); 
set_XY_Tick;

% subplot(2,1,2);
axes(ha(2)); 
plot(k_PFHI_bm,imagw_PFHI_bm,'-','Color',pltc(1,:),'linewidth',3);
hold on;
plot(k_OFHI_bm,imagw_OFHI_bm,'-','Color',pltc(2,:),'linewidth',3);
%
plot(k_wi_th0_Astfalk2017Fig1,wi_th0_Astfalk2017Fig1,'k--','LineWidth',3);
plot(k_wi_th45_Astfalk2017Fig1,wi_th45_Astfalk2017Fig1,'k--','LineWidth',3);
%
% pl = line([0,3],[0,0]);
pl.Color = 'black';
pl.LineStyle = '--';
%
xlabel('$kd_{i}$','Interpreter','latex','FontSize',32);
ylabel('$\omega_i/\omega_{ci}$','Interpreter','latex','FontSize',32);


% xlim([0 3]);
ylim([0.0,0.09]);
set_XY_Tick;


% save figure
savefig('benchmark_Asftalk2017Fig1.fig');

set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],...
  'PaperSize',[screenposition(3:4)]);

% print(gcf,['benchmark_Asftalk2017Fig1','.eps'],'-depsc');
print(gcf,['benchmark_Asftalk2017Fig1','.pdf'],'-dpdf');


