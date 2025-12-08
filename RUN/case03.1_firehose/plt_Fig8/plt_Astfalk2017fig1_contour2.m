% plot Astfalk2017fig1_contour.m
% baiweiphys@gmail.com

clear;
close all;
clc;

%% 1. BM
foldername = 'pbk_data';
fname = ['./',foldername,'/Astfalk2017Fig1_OFHI3d_bm.mat'];
load(fname);
par_bm = par;
ww_bm = ww;
wcs_bm = wcs;
[THETA_bm,KK0_bm] = meshgrid(theta,kk0);
%
ww_contour_bm = ww_bm;
index1_bm = imag(ww_contour_bm/wcs_bm(2)) < 5.0e-3;
index2_bm = imag(ww_contour_bm/wcs_bm(2)) > 1.0e0;
index_bm = index1_bm | index2_bm;
ww_contour_bm(index_bm) = NaN + 1i*NaN;
wr_contour_bm = real(ww_contour_bm)/wcs(2);
wi_contour_bm = imag(ww_contour_bm)/wcs(2);
num_contour_bm = 1;

%% 2.PBK
foldername = 'pbk_data';
fname = ['./',foldername,'/Astfalk2017Fig1_OFHI3d_pbk210.mat'];
load(fname);
par_pbk = par;
ww_pbk = ww;
wcs_pbk = wcs;
[THETA_pbk,KK0_pbk] = meshgrid(theta,kk0);
%
ww_contour_pbk = ww_pbk;
index1_pbk = imag(ww_contour_pbk/wcs_pbk(2)) < 5.0e-3;
index2_pbk = imag(ww_contour_pbk/wcs_pbk(2)) > 1.0e0;
index_pbk = index1_pbk | index2_pbk;
ww_contour_pbk(index_pbk) = NaN + 1i*NaN;
wr_contour_pbk = real(ww_contour_pbk)/wcs_pbk(2);
wi_contour_pbk = imag(ww_contour_pbk)/wcs_pbk(2);
num_contour_pbk = 1;

%% Growth rate and frequency contour plot
legend_txt=["(a) BM($\kappa_{\parallel}=\kappa_{\perp}=2$)", ...
            "(b) PBK($\sigma=0.5$)"];

h=figure('unit','normalized','Position',[0.02 0.45 0.6 0.49],'DefaultAxesFontSize',23);
% h=figure('unit','normalized','Position',[0.01 0.45 0.7 0.65],'DefaultAxesFontSize',23);
ha = tight_subplot(1,2,[0.0 0.02],[0.16 0.09],[0.08 0.1]);
%
cmin = 5e-3;
cmax = 0.08;
axes(ha(1)); 
for ii=1:num_contour_bm
    % contourf(KK0,THETA*180/pi,imagw_contour(:,:,ii)/wcs(3));
    contourf(KK0_bm(:,:),THETA_bm(:,:)*180/pi,wi_contour_bm(:,:,ii),30);
    hold on;
end
pl = line([0,0.8],[45,45]);
pl.Color = 'black';
pl.LineStyle = '--';
pl.LineWidth = 1;
%
caxis([cmin, cmax]);
% shading interp
colormap('jet');
% colorbar;
xlabel('$kd_i$','Interpreter','latex','FontSize',30);
ylabel('$\theta$','Interpreter','latex','FontSize',30);
% xlim([1e-3 1e5]);
% ylim([0 0.7]);
title(['(a) $\omega_i/\omega_{ci}$ ', '(BM, $J_{e,i}=8$)'],...
       'Interpreter','latex','FontSize',26);
set_XY_Tick;
% legend(legend_txt(2), 'Interpreter','latex','FontSize',20);
% legend("boxoff");
%

axes(ha(2)); 
for ii=1:num_contour_pbk
    % contourf(KK0,THETA*180/pi,imagw_contour(:,:,ii)/wcs(3));
    contourf(KK0_pbk(:,:),THETA_pbk(:,:)*180/pi,wi_contour_pbk(:,:,ii),30);
    hold on;
end
pl = line([0,0.8],[45,45]);
pl.Color = 'black';
pl.LineStyle = '--';
pl.LineWidth = 1;
%
set(ha(2),'YTickLabel','');
caxis([cmin, cmax]);
% shading interp
colormap('jet');
% colorbar;
xlabel('$kd_i$','Interpreter','latex','FontSize',30);
% ylabel('\theta');
% xlim([1e-3 1e5]);
% ylim([0 0.7]);
title(['(b) $\omega_i/\omega_{ci}$ ', ...
       '(PBK, $\kappa_{\parallel i}=2$,  $\kappa_{\perp i}=10$)'],...
       'Interpreter','latex','FontSize',26);
set_XY_Tick;

% legend(legend_txt(2), 'Interpreter','latex','FontSize',20);
% legend("boxoff");
% colorbar('Position', [0.91, 0.25, 0.015, 0.65]);  % [左, 下, 宽, 高]
colorbar('Position', [0.91, 0.195, 0.015, 0.68]);  % [左, 下, 宽, 高]

%% save
set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],...
  'PaperSize',[screenposition(3:4)]);

filename = "Astfalk2017Fig1_contour";
print(gcf,filename,'-dpdf');

% print(filename,'-dpdf','-bestfit');
