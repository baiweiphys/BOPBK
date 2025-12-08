% plt_benchmark_Cattaert2007Fig1.m

close all;
clear
% clc;

load('./pbkData/saveData_Cattaert2007Fig1km1.mat');
bopbk_kk_km1 = pas;
for ii=1:4
    bopbk_ww_km1(:,ii) = wws(:,1,ii);
end
bopbk_wr_km1 = real(bopbk_ww_km1);
bopbk_wi_km1 = imag(bopbk_ww_km1);

%
% load('./pbkData/saveData_Cattaert2007Fig1pbk1.mat');
% bopbk_kk_pbk1 = pas;
% for ii=1:4
%     bopbk_ww_pbk1(:,ii) = wws(:,1,ii);
% end
% bopbk_wr_pbk1 = real(bopbk_ww_pbk1);
% bopbk_wi_pbk1 = imag(bopbk_ww_pbk1);

%%
% load csv file
wi_km1 = csvread('./Cattaert2007Fig1/wi_Cattaert2007Fig1.csv',1,0);
wr_km1 = csvread('./Cattaert2007Fig1/wr_Cattaert2007Fig1.csv',1,0);

wr_x = wr_km1(:,1);
wr_upperRX = wr_km1(:,2);
wr_OP = wr_km1(:,3);
wr_LX = wr_km1(:,4);
wr_lowerRX = wr_km1(:,5);


wi_x = wi_km1(:,1);
wi_upperRX = wi_km1(:,2);
wi_lowerRX = wi_km1(:,3);
wi_LX = wi_km1(:,4);
wi_OP = wi_km1(:,5);


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
% rootNo = 2:8;
% elem = 3;
h=figure('unit','normalized','Position',[0.01 0.45 0.40 0.63],'DefaultAxesFontSize',22);
ha = tight_subplot(2,1,[0.03 .03],[0.12 .02],[0.19 0.05]);
% subplot(2,1,1);
axes(ha(1)); 
plot(bopbk_kk_km1,bopbk_wr_km1(:,1),'-','Color',pltc(1,:),'linewidth',3); % upper R-X
hold on;
plot(bopbk_kk_km1,bopbk_wr_km1(:,2),'-','Color',pltc(2,:),'linewidth',3); % O-P
plot(bopbk_kk_km1,bopbk_wr_km1(:,3),'-','Color',pltc(3,:),'linewidth',3); % L-X
plot(bopbk_kk_km1,bopbk_wr_km1(:,4),'-','Color',pltc(4,:),'linewidth',3); % lower R- X
%
plot(wr_x,wr_upperRX,'k--','Markersize',10,'linewidth',3);
plot(wr_x,wr_OP,'k--','Markersize',10,'linewidth',3);
plot(wr_x,wr_LX,'k--','Markersize',10,'linewidth',3);
plot(wr_x,wr_lowerRX,'k--','MarkerSize',10,'linewidth',3);
%
% plot(bopbk_kk_pbk1,bopbk_wr_pbk1(:,1),'--','Color',pltc(1,:),'linewidth',3); % upper R-X
% plot(bopbk_kk_pbk1,bopbk_wr_pbk1(:,2),'--','Color',pltc(2,:),'linewidth',3); % O-P
% plot(bopbk_kk_pbk1,bopbk_wr_pbk1(:,3),'--','Color',pltc(3,:),'linewidth',3); % L-X
% plot(bopbk_kk_pbk1,bopbk_wr_pbk1(:,4),'--','Color',pltc(4,:),'linewidth',3); % lower R-X
ylabel('$\omega_r/|\omega_{ce}|$','Interpreter','latex','FontSize',32);
% xlabel('k\rho_{ce}');
% xlim([1e-3 1e5]);
ylim([0 3]);
% legend('Upper R-X branch','O-P wave','L-X wave','Lower R-X branch','Cattaert2007',...
%        'Interpreter','latex','FontSize',20,'Location','northwest'); 
legend('Upper R-X branch','O-P wave','L-X wave','Lower R-X branch', ...
        'Cattaert(2007)',...
       'Interpreter','latex','FontSize',20,'Location','northwest'); 
legend('boxoff');
set(ha(1),'XTickLabel',''); 
set_XY_Tick;

%
% subplot(2,1,2);
axes(ha(2)); 
plot(bopbk_kk_km1,bopbk_wi_km1(:,1),'-','Color',pltc(1,:),'linewidth',3); % upper R-X
hold on;
plot(bopbk_kk_km1,bopbk_wi_km1(:,2),'-','Color',pltc(2,:),'linewidth',3); % O-P
plot(bopbk_kk_km1,bopbk_wi_km1(:,4),'-','Color',pltc(4,:),'linewidth',3); % L-X
plot(bopbk_kk_km1,bopbk_wi_km1(:,3),'-','Color',pltc(3,:),'linewidth',3); % lower R- X
%
plot(wi_x,wi_upperRX,'k--','Markersize',10,'linewidth',3);
plot(wi_x,wi_OP,'k--','Markersize',10,'linewidth',3);
plot(wi_x,wi_LX,'k--','Markersize',10,'linewidth',3);
plot(wi_x,wi_lowerRX,'k--','MarkerSize',10,'linewidth',3);
%
% plot(bopbk_kk_pbk1,bopbk_wi_pbk1(:,1),'--','Color',pltc(1,:),'linewidth',3); % upper R-X
% plot(bopbk_kk_pbk1,bopbk_wi_pbk1(:,2),'--','Color',pltc(2,:),'linewidth',3); % O-P
% plot(bopbk_kk_pbk1,bopbk_wi_pbk1(:,3),'--','Color',pltc(3,:),'linewidth',3); % L-X
% plot(bopbk_kk_pbk1,bopbk_wi_pbk1(:,4),'--','Color',pltc(4,:),'linewidth',3); % lower R-X
ylabel('$\omega_i/|\omega_{ce}|$','Interpreter','latex','FontSize',32);
xlabel('$k\rho_{ce}$','Interpreter','latex','FontSize',32);
% xlim([0,0.37]);
ylim([-0.01,0.0]);
set_XY_Tick;


% save figure
savefig('benchmark_Cattaert2007Fig1.fig');

set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],...
  'PaperSize',[screenposition(3:4)]);

% print(gcf,['benchmark_Cattaert2007Fig1','.eps'],'-depsc');
print(gcf,['benchmark_Cattaert2007Fig1','.pdf'],'-dpdf');
