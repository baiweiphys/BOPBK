% @Description: Benchmark cases visualization and saving.
% @Filename: plotall.m
% @Author: Bai Wei (baiwei12@mail.ustc.edu.cn, baiweiphys@gmail.com)
% Created: 2026-01-27
% Last Modified: 2026-01-27

clear all;
clc;
close all;

%% ================= LOG INITIALIZATION =================
addpath('./tools');
scriptName = 'plotall';  % Name of current script
logFile = initLogFile(scriptName);
diary(logFile);

% Log startup information
fprintf('Log file: %s\n', logFile);
fprintf('Start time: %s\n\n', datetime('now'));

%% ==================== MAIN PROGRAM ====================
disp('===========================================================');
disp('BO-PBK SOLVER');
disp('Dispersion Relation for Oblique Waves in Magnetized Plasma');
disp('with Anisotropic Drift Loss-Cone Distributions');
disp('===========================================================');
fprintf('Version: 1.0 beta\n');
fprintf('Date: %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
disp('===========================================================');
fprintf('\n');

% Plot initialization
fprintf('Initializing benchmark visualization...\n');


%% plt_Fig3
clear; close all;
%
plt_path = './results/plt_Fig3';
pltFile = 'plt_benchmark_Cattaert2007Fig1.m';
%
plot_script3 = fullfile(plt_path, pltFile);
fprintf('Running: %s\n', plot_script3);
%
if exist(plot_script3, 'file')
    run(plot_script3);
else
    warning('Script not found: %s', plot_script3);
end

%% plt_Fig4
clear; close all;
%
plt_path = './results/plt_Fig4';
pltFile = 'plt_benchmark_lazar2010Fig2.m';
%
plot_script4 = fullfile(plt_path, pltFile);
fprintf('Running: %s\n', plot_script4);
%
if exist(plot_script4, 'file')
    run(plot_script4);
else
    warning('Script not found: %s', plot_script4);
end


%% plt_Fig5
clear; close all;
%
plt_path = './results/plt_Fig5';
pltFile = 'plt_benchmark_Cattaert2007Fig9.m';
%
plot_script5 = fullfile(plt_path, pltFile);
fprintf('Running: %s\n', plot_script5);
%
if exist(plot_script5, 'file')
    run(plot_script5);
else
    warning('Script not found: %s', plot_script5);
end


%% plt_Fig6
clear; close all;
%
plt_path = './results/plt_Fig6';
pltFile = 'plt_benchmark_Astfalk2017Fig1.m';
%
plot_script6 = fullfile(plt_path, pltFile);
fprintf('Running: %s\n', plot_script6);
%
if exist(plot_script6, 'file')
    run(plot_script6);
else
    warning('Script not found: %s', plot_script6);
end

%% plt_Fig7
clear; close all;
%
plt_path = './results/plt_Fig7';
pltFile = 'plt_OFHI.m';
%
plot_script7 = fullfile(plt_path, pltFile);
fprintf('Running: %s\n', plot_script7);
%
if exist(plot_script7, 'file')
    run(plot_script7);
else
    warning('Script not found: %s', plot_script7);
end

%% plt_Fig8

% clear; close all;
% %
% plt_path = './results/plt_Fig8';
% pltFile = 'plt_Astfalk2017fig1_contour.m';
% %
% plot_script = fullfile(plt_path, pltFile);
% fprintf('Running: %s\n', plot_script);
% %
% if exist(plot_script, 'file')
%     run(plot_script);
% else
%     warning('Script not found: %s', plot_script);
% end


%% plt_Fig9
clear; close all;
%
plt_path = './results/plt_Fig9';
pltFile = 'plt_benchmark_Santos2014PoPFig6.m';
%
plot_script9 = fullfile(plt_path, pltFile);
fprintf('Running: %s\n', plot_script9);
%
if exist(plot_script9, 'file')
    run(plot_script9);
else
    warning('Script not found: %s', plot_script9);
end
%% plt_Fig10-11
clear; close all;
%
plt_path = './results/plt_Fig10-11';
pltFile10 = 'plt_Fig10.m';
pltFile11 = 'plt_Fig11.m';
%
plot_script10 = fullfile(plt_path, pltFile10);
plot_script11 = fullfile(plt_path, pltFile11);
%
fprintf('Running: %s\n', plot_script10);
fprintf('Running: %s\n', plot_script11);
%
if exist(plot_script10, 'file') && exist(plot_script11, 'file')
    run(plot_script10);
    run(plot_script11);
else
    warning('Script not found: %s or %s', plot_script10, plot_script11);
end

%% plt_Fig12
clear; close all;
%
plt_path = './results/plt_Fig12';
pltFile = 'plt_Basu2011PoP_Fig7_with_sig.m';
%
plot_script12 = fullfile(plt_path, pltFile);
fprintf('Running: %s\n', plot_script12);
%
if exist(plot_script12, 'file')
    run(plot_script12);
else
    warning('Script not found: %s', plot_script12);
end

%% plt_Fig13
clear; close all;
%
plt_path = './results/plt_Fig13';
pltFile = 'plt_benchmark_Santos2015PoPFig6.m';
%
plot_script13 = fullfile(plt_path, pltFile);
fprintf('Running: %s\n', plot_script13);
%
if exist(plot_script13, 'file')
    run(plot_script13);
else
    warning('Script not found: %s', plot_script13);
end

%% plt_Fig14
clear; close all;
%
plt_path = './results/plt_Fig14';
pltFile = 'plt_EMIC.m';
%
plot_script14 = fullfile(plt_path, pltFile);
fprintf('Running: %s\n', plot_script14);
%
if exist(plot_script14, 'file')
    run(plot_script14);
else
    warning('Script not found: %s', plot_script14);
end

clear; close all;


%% ================== LOG COMPLETION ====================
fprintf('Done at: %s\n', datetime('now'));

% Stop logging
diary off;

rmpath('./tools');
