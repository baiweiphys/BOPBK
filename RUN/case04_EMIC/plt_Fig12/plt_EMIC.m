% plt_EMIC.m

close all;
clear;
% clc;


%% PBK data
load('./pbkData/pbk_EMIC_th=0_B0=100nT_Tsz=1879.mat');
for ii=1:4
    ww(:,ii) = wws(:,1,ii);
end
k_pbk_th0 = pas;
wr_pbk_th0 = real(ww);
wi_pbk_th0 = imag(ww);

%
load('./pbkData/pbk_EMIC_th=15_B0=100nT_Tsz=1879.mat');
for ii=1:4
    ww(:,ii) = wws(:,1,ii);
end
k_pbk_th15 = pas;
wr_pbk_th15 = real(ww);
wi_pbk_th15 = imag(ww);

%
load('./pbkData/pbk_EMIC_th=40_B0=100nT_Tsz=1879.mat');
for ii=1:4
    ww(:,ii) = wws(:,1,ii);
end
k_pbk_th40 = pas;
wr_pbk_th40 = real(ww);
wi_pbk_th40 = imag(ww);

load('./pbkData/pbk_EMIC_th=55_B0=100nT_Tsz=1879.mat');
for ii=1:4
    ww(:,ii) = wws(:,1,ii);
end
k_pbk_th55 = pas;
wr_pbk_th55 = real(ww);
wi_pbk_th55 = imag(ww);

%% BOKM data
load('./bai2025CPCfig13_bokmData/EMIC_th=0_B0=100nT_Tsz=1879.mat');
for ii=1:4
    ww(:,ii) = wws(:,1,ii);
end
k_bokm_th0 = pas;
wr_bokm_th0 = real(ww);
wi_bokm_th0 = imag(ww);

%
load('./bai2025CPCfig13_bokmData/EMIC_th=15_B0=100nT_Tsz=1879.mat');
for ii=1:4
    ww(:,ii) = wws(:,1,ii);
end
k_bokm_th15 = pas;
wr_bokm_th15 = real(ww);
wi_bokm_th15 = imag(ww);

%
load('./bai2025CPCfig13_bokmData/EMIC_th=40_B0=100nT_Tsz=1879.mat');
for ii=1:4
    ww(:,ii) = wws(:,1,ii);
end
k_bokm_th40 = pas;
wr_bokm_th40 = real(ww);
wi_bokm_th40 = imag(ww);

%
load('./bai2025CPCfig13_bokmData/EMIC_th=55_B0=100nT_Tsz=1879.mat');
for ii=1:4
    ww(:,ii) = wws(:,1,ii);
end
k_bokm_th55 = pas;
wr_bokm_th55 = real(ww);
wi_bokm_th55 = imag(ww);

%%
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
h=figure('unit','normalized','Position',[0.01 0.45 0.8 0.55],'DefaultAxesFontSize',23);
ha = tight_subplot(2,4,[0.02 0.02],[0.15 0.07],[0.08 0.03]);

% fig(a)
axes(ha(1)); 
for ii=1:4
    plot(k_pbk_th0,wr_pbk_th0(:,ii),'-','Color',pltc(ii,:),'linewidth',3); %
    hold on;
end
%
for ii=1:4
    plot(k_bokm_th0,wr_bokm_th0(:,ii),'k--','linewidth',3); %
    hold on;
end
ylabel('$\omega_r/\omega_{cp}$','Interpreter','latex','FontSize',28);
xlim([0,0.08]);
ylim([0,1]);
pl = line([0,0.08],[0.456,0.456]);
pl.Color = 'black';
pl.LineStyle = '--';
%
pl = line([0,0.08],[0.185,0.185]);
pl.Color = 'black';
pl.LineStyle = '--';
%
legend('L-mode (O$^+$)','L-mode (He$^+$)','L-mode (H$^+$)','R-mode',...
       'BO-KM',...
       'Interpreter','latex','FontSize',15,'Location','northwest');
legend('boxoff');
title('(a) $\theta=0^{\circ}$','Interpreter','latex');
%
set_XY_Tick;
text(0.06,0.13,'$\omega_{cr1}/\omega_{cp}$','fontsize',15,'Interpreter','latex');
text(0.06,0.385,'$\omega_{cr2}/\omega_{cp}$','FontSize',15,'Interpreter','latex');

axes(ha(5)); 
% h=figure('unit','normalized','Position',[0.01 0.45 0.7 0.45],'DefaultAxesFontSize',20);
for ii=1:4
    plot(k_pbk_th0,wi_pbk_th0(:,ii),'-','Color',pltc(ii,:),'linewidth',3); %
    hold on;
end
%
for ii=1:4
    plot(k_bokm_th0,wi_bokm_th0(:,ii),'k--','linewidth',3); %
    hold on;
end
pl = line([0,0.08],[0,0]);
pl.Color = 'black';
pl.LineStyle = '--';
%
ylabel('$\omega_i/\omega_{cp}$','Interpreter','latex','FontSize',28);
xlabel('$k\rho_{cp}$','Interpreter','latex','FontSize',28);
xlim([0,0.08]);
ylim([-0.1,0.06]);
set_XY_Tick;

%%
% fig(b)
axes(ha(2)); 
for ii=1:4
    plot(k_pbk_th15,wr_pbk_th15(:,ii),'-','Color',pltc(ii,:),'linewidth',3); %
    hold on;
end
%
for ii=1:4
    plot(k_bokm_th15,wr_bokm_th15(:,ii),'k--','linewidth',3); %
    hold on;
end
xlim([0,0.08]);
ylim([0,1]);
pl = line([0,0.08],[0.456,0.456]);
pl.Color = 'black';
pl.LineStyle = '--';
%
pl = line([0,0.08],[0.185,0.185]);
pl.Color = 'black';
pl.LineStyle = '--';
%
legend('L-mode (O$^+$)','L-mode (He$^+$)','L-mode (H$^+$)','R-mode',...
       'BO-KM',...
       'Interpreter','latex','FontSize',15,'Location','northwest');
legend('boxoff');
title('(b) $\theta=15^{\circ}$','Interpreter','latex');
set_XY_Tick;

%
axes(ha(6)); 
% h=figure('unit','normalized','Position',[0.01 0.45 0.7 0.45],'DefaultAxesFontSize',20);
for ii=1:4
    plot(k_pbk_th15,wi_pbk_th15(:,ii),'-','Color',pltc(ii,:),'linewidth',3); %
    hold on;
end
%
for ii=1:4
    plot(k_bokm_th15,wi_bokm_th15(:,ii),'k--','linewidth',3); %
    hold on;
end
pl = line([0,0.08],[0,0]);
pl.Color = 'black';
pl.LineStyle = '--';
%
%ylabel('\gamma/|\omega_{cp}|');
xlabel('$k\rho_{cp}$','Interpreter','latex','FontSize',28);
xlim([0,0.08]);
ylim([-0.1,0.06]);
set_XY_Tick;

%%
% fig(c)
axes(ha(3)); 
for ii=1:4
    plot(k_pbk_th40,wr_pbk_th40(:,ii),'-','Color',pltc(ii,:),'linewidth',3); %
    hold on;
end
%
for ii=1:4
    plot(k_bokm_th40,wr_bokm_th40(:,ii),'k--','linewidth',3); %
    hold on;
end
xlim([0,0.08]);
ylim([0,1]);
%
pl = line([0,0.08],[0.456,0.456]);
pl.Color = 'black';
pl.LineStyle = '--';
%
pl = line([0,0.08],[0.185,0.185]);
pl.Color = 'black';
pl.LineStyle = '--';
%
legend('L-mode (O$^+$)','Branch 1','Branch 2','Branch 3',...
       'BO-KM',...
       'Interpreter','latex','FontSize',15,'Location','northwest');
legend('boxoff');
%
title('(c) $\theta=40^{\circ}$','Interpreter','latex');
set_XY_Tick;

%
axes(ha(7)); 
% h=figure('unit','normalized','Position',[0.01 0.45 0.7 0.45],'DefaultAxesFontSize',20);
for ii=1:4
    plot(k_pbk_th40,wi_pbk_th40(:,ii),'-','Color',pltc(ii,:),'linewidth',3); %
    hold on;
end
%
for ii=1:4
    plot(k_bokm_th40,wi_bokm_th40(:,ii),'k--','linewidth',3); %
    hold on;
end
pl = line([0,0.08],[0,0]);
pl.Color = 'black';
pl.LineStyle = '--';
%
xlabel('$k\rho_{cp}$','Interpreter','latex','FontSize',28);
xlim([0,0.08]);
ylim([-0.1,0.06]);
set_XY_Tick;


%%
% fig(d)
axes(ha(4)); 
for ii=1:4
    plot(k_pbk_th55,wr_pbk_th55(:,ii),'-','Color',pltc(ii,:),'linewidth',3); %
    hold on;
end
%
for ii=1:4
    plot(k_bokm_th55,wr_bokm_th55(:,ii),'k--','linewidth',3); %
    hold on;
end
%
pl = line([0,0.08],[0.456,0.456]);
pl.Color = 'black';
pl.LineStyle = '--';
%
pl = line([0,0.08],[0.185,0.185]);
pl.Color = 'black';
pl.LineStyle = '--';
%
legend('L-mode (O$^+$)','Branch 1','Branch 2','Branch 3',...
       'BO-KM',...
       'Interpreter','latex','FontSize',15,'Location','northwest');
legend('boxoff');
xlim([0,0.08]);
ylim([0,1]);
title('(d) $\theta=55^{\circ}$','Interpreter','latex');
set_XY_Tick;

%
axes(ha(8)); 
% h=figure('unit','normalized','Position',[0.01 0.45 0.7 0.45],'DefaultAxesFontSize',20);
for ii=1:4
    plot(k_pbk_th55,wi_pbk_th55(:,ii),'-','Color',pltc(ii,:),'linewidth',3); %
    hold on;
end
%
for ii=1:4
    plot(k_bokm_th55,wi_bokm_th55(:,ii),'k--','linewidth',3); %
    hold on;
end
pl = line([0,0.08],[0,0]);
pl.Color = 'black';
pl.LineStyle = '--';
%
xlabel('$k\rho_{cp}$','Interpreter','latex','FontSize',28);
xlim([0,0.08]);
ylim([-0.1,0.06]);
set_XY_Tick;

%%
for ii=1:8
    if(ii==1 | ii==5)
        ;
    else
        set(ha(ii),'YTickLabel','');
    end

    if (ii<5)
        set(ha(ii),'XTickLabel',''); 
    end
end


%% save figure
savefig('EMIC_B0=100nT_Tsz=1879.fig');

set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],...
  'PaperSize',[screenposition(3:4)]);

% print(gcf,['EMIC_B0=100nT_Tsz=1879','.eps'],'-depsc');
print(gcf,['EMIC_B0=100nT_Tsz=1879','.pdf'],'-dpdf');
