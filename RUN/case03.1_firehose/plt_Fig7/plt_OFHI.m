
clear;
clc;

%% read bopbk data
load('./pbk_data/saveData_Astfalk2017Fig1_OFHI_bm.mat');
k_OFHI_bm = pas;
w_OFHI_bm = wws;
realw_OFHI_bm = real(w_OFHI_bm);
imagw_OFHI_bm = imag(w_OFHI_bm);

load('./pbk_data/saveData_Astfalk2017Fig1_OFHI_pbk2.mat');
k_OFHI_pbk2 = pas;
w_OFHI_pbk2 = wws;
realw_OFHI_pbk2 = real(w_OFHI_pbk2);
imagw_OFHI_pbk2 = imag(w_OFHI_pbk2);

load('./pbk_data/saveData_Astfalk2017Fig1_OFHI_pbk4.mat');
k_OFHI_pbk4 = pas;
w_OFHI_pbk4 = wws;
realw_OFHI_pbk4 = real(w_OFHI_pbk4);
imagw_OFHI_pbk4 = imag(w_OFHI_pbk4);

load('./pbk_data/saveData_Astfalk2017Fig1_OFHI_pbk8.mat');
k_OFHI_pbk8 = pas;
w_OFHI_pbk8 = wws;
realw_OFHI_pbk8 = real(w_OFHI_pbk8);
imagw_OFHI_pbk8 = imag(w_OFHI_pbk8);

%
load('./pbk_data/saveData_Astfalk2017Fig1_OFHI_km2.mat');
k_OFHI_km2 = pas;
w_OFHI_km2 = wws;
realw_OFHI_km2 = real(w_OFHI_km2);
imagw_OFHI_km2 = imag(w_OFHI_km2);

load('./pbk_data/saveData_Astfalk2017Fig1_OFHI_km4.mat');
k_OFHI_km4 = pas;
w_OFHI_km4 = wws;
realw_OFHI_km4 = real(w_OFHI_km4);
imagw_OFHI_km4 = imag(w_OFHI_km4);

load('./pbk_data/saveData_Astfalk2017Fig1_OFHI_km8.mat');
k_OFHI_km8 = pas;
w_OFHI_km8 = wws;
realw_OFHI_km8 = real(w_OFHI_km8);
imagw_OFHI_km8 = imag(w_OFHI_km8);

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

%% plot 
% h=figure('unit','normalized','Position',[0.01 0.45 0.4 0.54],'DefaultAxesFontSize',14);
h=figure('unit','normalized','Position',[0.01 0.45 0.42 0.45],'DefaultAxesFontSize',22);
ha = tight_subplot(1,1,[0.03 .03],[0.18 0.03],[0.14 0.03]);
axes(ha);
% PBK data
plot(k_OFHI_pbk2,imagw_OFHI_pbk2(:,1,1),'-','Color',pltc(1,:),'linewidth',3);
hold on;
% plot(k_OFHI_pbk2,imagw_OFHI_pbk2(:,1,2),'-','Color',pltc(1,:),'linewidth',3);
plot(k_OFHI_pbk4,imagw_OFHI_pbk4,'-','Color',pltc(2,:),'linewidth',3);
plot(k_OFHI_pbk8,imagw_OFHI_pbk8,'-','Color',pltc(3,:),'linewidth',3);
% KM data
plot(k_OFHI_km2,imagw_OFHI_km2,'-.','Color',pltc(1,:),'linewidth',3);
plot(k_OFHI_km4,imagw_OFHI_km4,'-.','Color',pltc(2,:),'linewidth',3);
plot(k_OFHI_km8,imagw_OFHI_km8,'-.','Color',pltc(3,:),'linewidth',3);
plot(k_OFHI_bm,imagw_OFHI_bm,'-','Color',pltc(8,:),'linewidth',3);
% grid on;
xlabel('$kd_{i}$','Interpreter','latex','FontSize',33);
ylabel('$\omega_i/\omega_{ci}$','Interpreter','latex','FontSize',33);
%
xlim([0,0.5])
ylim([0,0.08]);
legend('PBK($\kappa_{\parallel i}=2$,$\kappa_{\perp i}=10$)',...
       'PBK($\kappa_{\parallel i}=4$,$\kappa_{\perp i}=10$)',...
       'PBK($\kappa_{\parallel i}=8$,$\kappa_{\perp i}=10$)',...
       'KM($\kappa_{\parallel i}=2$)',...
       'KM($\kappa_{\parallel i}=4$)',...
       'KM($\kappa_{\parallel i}=8$)',...
       'BM($J_{e,i}=8$)',...
       'Interpreter','latex','FontSize',22,'Location','northwest');
legend('boxoff');

set_XY_Tick;

% save figure
savefig('benchmark_OFHI.fig');

set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],...
  'PaperSize',[screenposition(3:4)]);

% print(gcf,['benchmark_OFHI','.pdf'],'-depsc');
print(gcf,['benchmark_OFHI','.pdf'],'-dpdf');
