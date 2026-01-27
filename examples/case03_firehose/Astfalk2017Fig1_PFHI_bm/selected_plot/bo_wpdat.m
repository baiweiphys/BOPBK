% 18-10-19 17:56 Hua-sheng XIE, huashengxie@gmail.com, CCF-ENN, China
% Initial data for run bo_plot_select

% Search the most close dispersion surfaces to these data.
% Initial data for find the corresponding dispersion surfaces.
% Please use bo_plot_all.m to visualize all the solutions, and then
% modify here the initial point of which mode(s) you want plot/store.

% wpdat(:,1) is pa; wpdat(:,2) is pb for 2D scan and arbitrary for 1D scan;
% wpdat(:,3) is Re or Im(omega)


wpdat=[0.166103,0,0.210987i; 
       %0.0835517,0,1.11005;
       %0.152345,0,0.85164;
       %0.152345,0,0.694488;
       %0.138586,0,0.506415;
  ];


jselc=1; % jselc=1, alway plot the most unstable dispersion surface
% wpdat=[];

%%
% new add

npa = 100;
npb = 1;
nw = 849;
ipa = 1;
ipb = 1;
iloga = 0;
jsetplot = 1;
iem = 3;

wwn = www/wcs(2);
pb = 0;
pa = kk0;


%%
run ./bo_plot_select;
% subplot(122);xlim([0,2]);ylim([-1e-3,2e-3]);subplot(121);xlim([0,2]);
% subplot(122);ylim([-0.5,0.5]);
