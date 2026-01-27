% 18-10-19 17:56 Hua-sheng XIE, huashengxie@gmail.com, CCF-ENN, China
% Initial data for run bo_plot_select

% Search the most close dispersion surfaces to these data.
% Initial data for find the corresponding dispersion surfaces.
% Please use bo_plot_all.m to visualize all the solutions, and then
% modify here the initial point of which mode(s) you want plot/store.

% wpdat(:,1) is pa; wpdat(:,2) is pb for 2D scan and arbitrary for 1D scan;
% wpdat(:,3) is Re or Im(omega)

wpdat=[
       % 55, 0.0080898, 0.0339003;    % L mode O+ band
       % 55, 0.00889778, 0.0951292; % L mode He+ band
       % 55, 0.0339452, 0.326385;     % L mode proton band
       % 55, 0.0525287, 0.633827;   % R mode
       %
       40, 0.00647384, 0.0361281; % L mode O+ band
       40, 0.0315212, 0.384354;   % L mode proton band
       40, 0.0210175, 0.417502;   % L mode proton band  
       40, 0.0258654, 0.201974;     % R mode 
  ];
jselc=0; % alway plot the most unstable dispersion surface
% jselc = 0; 
% wpdat=[];

%%
% new add
wwn = ww/abs(wcs(1));
pa = theta*180/pi;
pb = kk*rhocs(1);
npa = size(ww,1);
npb = size(ww,2);
nw = size(ww,3);
ipa = 1;
ipb = 2;
iem = 3;

%%
run ./bo_plot_select;
% subplot(122);xlim([0,2]);ylim([-1e-3,2e-3]);subplot(121);xlim([0,2]);
% subplot(122);ylim([-0.5,0.5]);
