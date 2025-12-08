% 18-10-19 17:56 Hua-sheng XIE, huashengxie@gmail.com, CCF-ENN, China
% Initial data for run bo_plot_select

% Search the most close dispersion surfaces to these data.
% Initial data for find the corresponding dispersion surfaces.
% Please use bo_plot_all.m to visualize all the solutions, and then
% modify here the initial point of which mode(s) you want plot/store.

% wpdat(:,1) is pa; wpdat(:,2) is pb for 2D scan and arbitrary for 1D scan;
% wpdat(:,3) is Re or Im(omega)


wpdat=[ 
        % 0.0080898,  0, 0.0339003; % L mode O+ band
        % 0.00889778, 0, 0.0951292; % L mode He+ band
        % 0.0339452,  0, 0.326385; % L mode proton band
        % 0.0525287,  0, 0.633827; % R mode
        % 0.0153616,  0, 0.040865; 
  ];

% jselc=1, alway plot the most unstable dispersion surface
jselc=1; 
% wpdat=[];

%%
% new add

npa = 100;
npb = 1;
nw = 792;
ipa = 1;
ipb = 1;
iloga = 0;
jsetplot = 1;
iem = 3;

wwn = www/abs(wcs(1));
pb = 0;
pa = kk0;


%%
run ./bo_plot_select;
% subplot(122);xlim([0,2]);ylim([-1e-3,2e-3]);subplot(121);xlim([0,2]);
% subplot(122);ylim([-0.5,0.5]);
