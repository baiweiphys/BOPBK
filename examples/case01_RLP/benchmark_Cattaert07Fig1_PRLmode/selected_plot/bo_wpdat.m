% 18-10-19 17:56 Hua-sheng XIE, huashengxie@gmail.com, CCF-ENN, China
% Initial data for run bo_plot_select

% Search the most close dispersion surfaces to these data.
% Initial data for find the corresponding dispersion surfaces.
% Please use bo_plot_all.m to visualize all the solutions, and then
% modify here the initial point of which mode(s) you want plot/store.

% wpdat(:,1) is pa; wpdat(:,2) is pb for 2D scan and arbitrary for 1D scan;
% wpdat(:,3) is Re or Im(omega)


wpdat=[0.1,0,1.37366;
       0.1,0,1.07271;
       0.1,0,0.75519;
       0.1,0,0.39861;
       % 0.094,0,-0.00026i;
       % 0.106,0,-0.00104i;
       % 0.090979,0,-0.002076i;
       % 0.181858,0,-0.00475962i;
  ];


jselc=0; % jselc=1, alway plot the most unstable dispersion surface
% wpdat=[];

%%
% new add

npa = 200;
npb = 1;
nw = 114;
ipa = 1;
ipb = 1;
iloga = 0;
jsetplot = 1;
iem = 3;

wwn = www/abs(wcs(1));
pb = 0;
pa = kk*rhocs(1);


%%
run ./bo_plot_select;
% subplot(122);xlim([0,2]);ylim([-1e-3,2e-3]);subplot(121);xlim([0,2]);
% subplot(122);ylim([-0.5,0.5]);
