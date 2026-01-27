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
plot_script9 = fullfile(plt_path, pltFile);
fprintf('Running: %s\n', plot_script9);
%
if exist(plot_script9, 'file')
    run(plot_script9);
else
    warning('Script not found: %s', plot_script9);
end

%% plt_Fig4
clear; close all;
%
plt_path = './results/plt_Fig4';
pltFile = 'plt_benchmark_lazar2011Fig2.m';
%
plot_script9 = fullfile(plt_path, pltFile);
fprintf('Running: %s\n', plot_script9);
%
if exist(plot_script9, 'file')
    run(plot_script9);
else
    warning('Script not found: %s', plot_script9);
end


%% plt_Fig5
clear; close all;
%
plt_path = './results/plt_Fig5';
pltFile = 'plt_benchmark_Cattaert2007Fig9.m';
%
plot_script9 = fullfile(plt_path, pltFile);
fprintf('Running: %s\n', plot_script9);
%
if exist(plot_script9, 'file')
    run(plot_script9);
else
    warning('Script not found: %s', plot_script9);
end


%% plt_Fig6
clear; close all;
%
plt_path = './results/plt_Fig6';
pltFile = 'plt_benchmark_Astfalk2017Fig1.m';
%
plot_script9 = fullfile(plt_path, pltFile);
fprintf('Running: %s\n', plot_script9);
%
if exist(plot_script9, 'file')
    run(plot_script9);
else
    warning('Script not found: %s', plot_script9);
end

%% plt_Fig7
clear; close all;
%
plt_path = './results/plt_Fig7';
pltFile = 'plt_OFHI.m';
%
plot_script9 = fullfile(plt_path, pltFile);
fprintf('Running: %s\n', plot_script9);
%
if exist(plot_script9, 'file')
    run(plot_script9);
else
    warning('Script not found: %s', plot_script9);
end

%% plt_Fig8

% clear; close all;
% %
% plt_path = './results/plt_Fig8';
% pltFile = 'plt_Astfalk2017fig1_contour2.m';
% %
% plot_script = fullfile(plt_path, pltFile);
% fprintf('Running: %s\n', plot_script);
% %
% if exist(plot_script, 'file')
%     run(plot_script);
% else
%     warning('Script not found: %s', plot_script);
% end


%% plt_Fig9-10
clear; close all;
%
plt_path = './results/plt_Fig9-10';
pltFile9 = 'plt_Fig9.m';
pltFile10 = 'plt_Fig10.m';
%
plot_script9 = fullfile(plt_path, pltFile9);
plot_script10 = fullfile(plt_path, pltFile10);
%
fprintf('Running: %s\n', plot_script9);
fprintf('Running: %s\n', plot_script10);
%
if exist(plot_script9, 'file') && exist(plot_script10, 'file')
    run(plot_script9);
    run(plot_script10);
else
    warning('Script not found: %s or %s', plot_script9, plot_script10);
end


%% plt_Fig11
clear; close all;
%
plt_path = './results/plt_Fig11';
pltFile = 'plt_Basu2011PoP_Fig7_with_sig.m';
%
plot_script9 = fullfile(plt_path, pltFile);
fprintf('Running: %s\n', plot_script9);
%
if exist(plot_script9, 'file')
    run(plot_script9);
else
    warning('Script not found: %s', plot_script9);
end


%% plt_Fig12
clear; close all;
%
plt_path = './results/plt_Fig12';
pltFile = 'plt_EMIC.m';
%
plot_script9 = fullfile(plt_path, pltFile);
fprintf('Running: %s\n', plot_script9);
%
if exist(plot_script9, 'file')
    run(plot_script9);
else
    warning('Script not found: %s', plot_script9);
end

clear; close all;


%% ================== LOG COMPLETION ====================
fprintf('Done at: %s\n', datetime('now'));

% Stop logging
diary off;

rmpath('./tools');
