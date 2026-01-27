% @Description: To execute all benchmark cases in the "./examples" fiels.
% BO-PBK (Version 1.0 beta, available at the URL: 
% https://github.com/baiweiphys/BOPBK/) is a new tool developed 
% by Wei Bai (TYUT) and Huasheng Xie (ENN). 
% Its primary objective is to efficiently calculate all solutions 
% for the dispersion relation of obliquely propagating waves in 
% magnetized hot plasmas with multiple species. 
% The code supports a wide range of multi-component velocity distributions, 
% including:
% Anisotropic drift loss-cone product-bi-kappa (PBK);
% Anisotropic drift loss-cone kappa-Maxwellian (KM);
% Anisotropic drift loss-cone bi-Maxwellian (BM);
% Hybrid combinations of these distributions.
% @Author: Bai Wei (baiwei12@mail.ustc.edu.cn, baiweiphys@gmail.com)
% @Date: 2026-01-26
% @LastEditors: Bai Wei
% @LastEditTime: 2026-01-27

clear;
clc;
close all;

%% ================= LOG INITIALIZATION =================
addpath('./tools');
scriptName = 'runall';  % Name of current script
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

fprintf('Starting benchmark calculations...\n');


%% Case#1
clear; close all;
%
case_path = './examples/case01_RLP/benchmark_Cattaert07Fig1_PRLmode';
dataFile = 'saveData_Cattaert2007Fig1km1.mat';
%
main_script = fullfile(case_path, 'main_bopbk.m');
plot_script = fullfile(case_path, 'selected_plot','bo_wpdat.m');
pbkdata = fullfile(case_path, 'selected_plot', dataFile);
%
fprintf('\n========================================\n');
fprintf('Running: %s\n', main_script);
%
if exist(main_script, 'file') && exist(plot_script, 'file')
    run(main_script);
    run(plot_script);
else
    warning('Script not found: %s or %s', main_script, plot_script);
end

% Save results
currentPath = pwd;
datapath = './results/plt_Fig3/pbkData';
run('./tools/createDateFile(currentPath,datapath)');
copyfile(pbkdata,datapath);

%% Case#2.1（BM）
clear; close all;
%
case_path = './examples/case02.1_whistler';
subcase = 'Whistler_Lazar11Fig2_bm'; 
dataFile = 'saveData_Electron_Whistler_bm8_Lazar2011Fig2.mat';
%
main_script = fullfile(case_path,subcase,'main_bopbk.m');
plot_script = fullfile(case_path,subcase,'selected_plot','bo_wpdat.m');
pbkdata = fullfile(case_path,subcase,'selected_plot', dataFile);
%
fprintf('\n========================================\n');
fprintf('Running: %s\n', main_script);
%
if exist(main_script, 'file') && exist(plot_script, 'file')
    run(main_script);
    run(plot_script);
else
    warning('Script not found: %s or %s', main_script, plot_script);
end

% Save results
currentPath = pwd;
datapath = './results/plt_Fig4/pbkData';
run('./tools/createDateFile(currentPath,datapath)');
copyfile(pbkdata,datapath);

%% Case#2.1 (PBK2)
clear; close all;
%
case_path = './examples/case02.1_whistler';
subcase = 'Whistler_Lazar11Fig2_pbk2'; 
dataFile = 'saveData_Electron_Whistler_pbk2_Lazar2011Fig2.mat';
%
main_script = fullfile(case_path,subcase,'main_bopbk.m');
plot_script = fullfile(case_path,subcase,'selected_plot','bo_wpdat.m');
pbkdata = fullfile(case_path,subcase,'selected_plot', dataFile);
%
fprintf('\n========================================\n');
fprintf('Running: %s\n', main_script);
%
if exist(main_script, 'file') && exist(plot_script, 'file')
    run(main_script);
    run(plot_script);
else
    warning('Script not found: %s or %s', main_script, plot_script);
end

% Save results
currentPath = pwd;
datapath = './results/plt_Fig4/pbkData';
run('./tools/createDateFile(currentPath,datapath)');
copyfile(pbkdata,datapath);


%% Case#2.1 (PBK6)
clear; close all;
%
case_path = './examples/case02.1_whistler';
subcase = 'Whistler_Lazar11Fig2_pbk6'; 
dataFile = 'saveData_Electron_Whistler_pbk6_Lazar2011Fig2.mat';
%
main_script = fullfile(case_path,subcase,'main_bopbk.m');
plot_script = fullfile(case_path,subcase,'selected_plot','bo_wpdat.m');
pbkdata = fullfile(case_path,subcase,'selected_plot', dataFile);
%
fprintf('\n========================================\n');
fprintf('Running: %s\n', main_script);
%
if exist(main_script, 'file') && exist(plot_script, 'file')
    run(main_script);
    run(plot_script);
else
    warning('Script not found: %s or %s', main_script, plot_script);
end

% Save results
currentPath = pwd;
datapath = './results/plt_Fig4/pbkData';
run('./tools/createDateFile(currentPath,datapath)');
copyfile(pbkdata,datapath);

%% Case#2.2 (Tex/Tez=0.5)
clear; close all;
%
case_path = './examples/case02.2_whistler';
subcase = 'benchmark_Cattaert07Fig9_Whistler/Whistler_TexTez=0.5_km1'; 
dataFile = 'saveData_Cattaert2007Fig9_TexTez0dot5.mat';
%
main_script = fullfile(case_path,subcase,'main_bopbk.m');
plot_script = fullfile(case_path,subcase,'selected_plot','bo_wpdat.m');
pbkdata = fullfile(case_path,subcase,'selected_plot', dataFile);
%
fprintf('\n========================================\n');
fprintf('Running: %s\n', main_script);
%
if exist(main_script, 'file') && exist(plot_script, 'file')
    run(main_script);
    run(plot_script);
else
    warning('Script not found: %s or %s', main_script, plot_script);
end

% Save results
currentPath = pwd;
datapath = './results/plt_Fig5/pbkData';
run('./tools/createDateFile(currentPath,datapath)');
copyfile(pbkdata,datapath);

%% Case#2.2 (Tex/Tez=1)
clear; close all;
%
case_path = './examples/case02.2_whistler';
subcase = 'benchmark_Cattaert07Fig9_Whistler/Whistler_TexTez=1_km1'; 
dataFile = 'saveData_Cattaert2007Fig9_TexTez1.mat';
%
main_script = fullfile(case_path,subcase,'main_bopbk.m');
plot_script = fullfile(case_path,subcase,'selected_plot','bo_wpdat.m');
pbkdata = fullfile(case_path,subcase,'selected_plot', dataFile);
%
fprintf('\n========================================\n');
fprintf('Running: %s\n', main_script);
%
if exist(main_script, 'file') && exist(plot_script, 'file')
    run(main_script);
    run(plot_script);
else
    warning('Script not found: %s or %s', main_script, plot_script);
end

% Save results
currentPath = pwd;
datapath = './results/plt_Fig5/pbkData';
run('./tools/createDateFile(currentPath,datapath)');
copyfile(pbkdata,datapath);

%% Case#2.2 (Tex/Tez=2)
clear; close all;
%
case_path = './examples/case02.2_whistler';
subcase = 'benchmark_Cattaert07Fig9_Whistler/Whistler_TexTez=2_km1'; 
dataFile = 'saveData_Cattaert2007Fig9_TexTez2.mat';
%
main_script = fullfile(case_path,subcase,'main_bopbk.m');
plot_script = fullfile(case_path,subcase,'selected_plot','bo_wpdat.m');
pbkdata = fullfile(case_path,subcase,'selected_plot', dataFile);
%
fprintf('\n========================================\n');
fprintf('Running: %s\n', main_script);
%
if exist(main_script, 'file') && exist(plot_script, 'file')
    run(main_script);
    run(plot_script);
else
    warning('Script not found: %s or %s', main_script, plot_script);
end

% Save results
currentPath = pwd;
datapath = './results/plt_Fig5/pbkData';
run('./tools/createDateFile(currentPath,datapath)');
copyfile(pbkdata,datapath);

%% Case#2.2 (Tex/Tez=3)
clear; close all;
%
case_path = './examples/case02.2_whistler';
subcase = 'benchmark_Cattaert07Fig9_Whistler/Whistler_TexTez=3_km1'; 
dataFile = 'saveData_Cattaert2007Fig9_TexTez3.mat';
%
main_script = fullfile(case_path,subcase,'main_bopbk.m');
plot_script = fullfile(case_path,subcase,'selected_plot','bo_wpdat.m');
pbkdata = fullfile(case_path,subcase,'selected_plot', dataFile);
%
fprintf('\n========================================\n');
fprintf('Running: %s\n', main_script);
%
if exist(main_script, 'file') && exist(plot_script, 'file')
    run(main_script);
    run(plot_script);
else
    warning('Script not found: %s or %s', main_script, plot_script);
end

% Save results
currentPath = pwd;
datapath = './results/plt_Fig5/pbkData';
run('./tools/createDateFile(currentPath,datapath)');
copyfile(pbkdata,datapath);

%% Case#3 (OFHI_bm)
clear; close all;
%
case_path = './examples/case03_firehose';
subcase = 'Astfalk2017Fig1_OFHI_bm'; 
dataFile = 'saveData_Astfalk2017Fig1_OFHI_bm.mat';
%
main_script = fullfile(case_path,subcase,'main_bopbk.m');
plot_script = fullfile(case_path,subcase,'selected_plot','bo_wpdat.m');
pbkdata = fullfile(case_path,subcase,'selected_plot', dataFile);
%
fprintf('\n========================================\n');
fprintf('Running: %s\n', main_script);
%
if exist(main_script, 'file') && exist(plot_script, 'file')
    run(main_script);
    run(plot_script);
else
    warning('Script not found: %s or %s', main_script, plot_script);
end

% Save results
currentPath = pwd;
datapath = './results/plt_Fig6/pbkData';
run('./tools/createDateFile(currentPath,datapath)');
copyfile(pbkdata,datapath);

% Save results to plt_Fig7
datapath = './results/plt_Fig7/pbkData';
run('./tools/createDateFile(currentPath,datapath)');
copyfile(pbkdata,datapath);

%% Case#3 (PFHI_bm)
clear; close all;
%
case_path = './examples/case03_firehose';
subcase = 'Astfalk2017Fig1_PFHI_bm'; 
dataFile = 'saveData_Astfalk2017Fig1_PFHI_bm.mat';
%
main_script = fullfile(case_path,subcase,'main_bopbk.m');
plot_script = fullfile(case_path,subcase,'selected_plot','bo_wpdat.m');
pbkdata = fullfile(case_path,subcase,'selected_plot', dataFile);
%
fprintf('\n========================================\n');
fprintf('Running: %s\n', main_script);
%
if exist(main_script, 'file') && exist(plot_script, 'file')
    run(main_script);
    run(plot_script);
else
    warning('Script not found: %s or %s', main_script, plot_script);
end

% Save results
currentPath = pwd;
datapath = './results/plt_Fig6/pbkData';
run('./tools/createDateFile(currentPath,datapath)');
copyfile(pbkdata,datapath);

%% Case#3 (OFHI_km2)
clear; close all;
%
case_path = './examples/case03_firehose';
subcase = 'Astfalk2017Fig1_OFHI_km2'; 
dataFile = 'saveData_Astfalk2017Fig1_OFHI_km2.mat';
%
main_script = fullfile(case_path,subcase,'main_bopbk.m');
plot_script = fullfile(case_path,subcase,'selected_plot','bo_wpdat.m');
pbkdata = fullfile(case_path,subcase,'selected_plot', dataFile);
%
fprintf('\n========================================\n');
fprintf('Running: %s\n', main_script);
%
if exist(main_script, 'file') && exist(plot_script, 'file')
    run(main_script);
    run(plot_script);
else
    warning('Script not found: %s or %s', main_script, plot_script);
end

% Save results
currentPath = pwd;
datapath = './results/plt_Fig7/pbkData';
run('./tools/createDateFile(currentPath,datapath)');
copyfile(pbkdata,datapath);

%% Case#3 (OFHI_km4)
clear; close all;
%
case_path = './examples/case03_firehose';
subcase = 'Astfalk2017Fig1_OFHI_km4'; 
dataFile = 'saveData_Astfalk2017Fig1_OFHI_km4.mat';
%
main_script = fullfile(case_path,subcase,'main_bopbk.m');
plot_script = fullfile(case_path,subcase,'selected_plot','bo_wpdat.m');
pbkdata = fullfile(case_path,subcase,'selected_plot', dataFile);
%
fprintf('\n========================================\n');
fprintf('Running: %s\n', main_script);
%
if exist(main_script, 'file') && exist(plot_script, 'file')
    run(main_script);
    run(plot_script);
else
    warning('Script not found: %s or %s', main_script, plot_script);
end

% Save results
currentPath = pwd;
datapath = './results/plt_Fig7/pbkData';
run('./tools/createDateFile(currentPath,datapath)');
copyfile(pbkdata,datapath);


%% Case#3 (OFHI_km8)
clear; close all;
%
case_path = './examples/case03_firehose';
subcase = 'Astfalk2017Fig1_OFHI_km8'; 
dataFile = 'saveData_Astfalk2017Fig1_OFHI_km8.mat';
%
main_script = fullfile(case_path,subcase,'main_bopbk.m');
plot_script = fullfile(case_path,subcase,'selected_plot','bo_wpdat.m');
pbkdata = fullfile(case_path,subcase,'selected_plot', dataFile);
%
fprintf('\n========================================\n');
fprintf('Running: %s\n', main_script);
%
if exist(main_script, 'file') && exist(plot_script, 'file')
    run(main_script);
    run(plot_script);
else
    warning('Script not found: %s or %s', main_script, plot_script);
end

% Save results
currentPath = pwd;
datapath = './results/plt_Fig7/pbkData';
run('./tools/createDateFile(currentPath,datapath)');
copyfile(pbkdata,datapath);


%% Case#3 (OFHI_pbk2)
clear; close all;
%
case_path = './examples/case03_firehose';
subcase = 'Astfalk2017Fig1_OFHI_pbk2'; 
dataFile = 'saveData_Astfalk2017Fig1_OFHI_pbk2.mat';
%
main_script = fullfile(case_path,subcase,'main_bopbk.m');
plot_script = fullfile(case_path,subcase,'selected_plot','bo_wpdat.m');
pbkdata = fullfile(case_path,subcase,'selected_plot', dataFile);
%
fprintf('\n========================================\n');
fprintf('Running: %s\n', main_script);
%
if exist(main_script, 'file') && exist(plot_script, 'file')
    run(main_script);
    run(plot_script);
else
    warning('Script not found: %s or %s', main_script, plot_script);
end

% Save results
currentPath = pwd;
datapath = './results/plt_Fig7/pbkData';
run('./tools/createDateFile(currentPath,datapath)');
copyfile(pbkdata,datapath);


%% Case#3 (OFHI_pbk4)
clear; close all;
%
case_path = './examples/case03_firehose';
subcase = 'Astfalk2017Fig1_OFHI_pbk4'; 
dataFile = 'saveData_Astfalk2017Fig1_OFHI_pbk4.mat';
%
main_script = fullfile(case_path,subcase,'main_bopbk.m');
plot_script = fullfile(case_path,subcase,'selected_plot','bo_wpdat.m');
pbkdata = fullfile(case_path,subcase,'selected_plot', dataFile);
%
fprintf('\n========================================\n');
fprintf('Running: %s\n', main_script);
%
if exist(main_script, 'file') && exist(plot_script, 'file')
    run(main_script);
    run(plot_script);
else
    warning('Script not found: %s or %s', main_script, plot_script);
end

% Save results
currentPath = pwd;
datapath = './results/plt_Fig7/pbkData';
run('./tools/createDateFile(currentPath,datapath)');
copyfile(pbkdata,datapath);


%% Case#3 (OFHI_pbk8)
clear; close all;
%
case_path = './examples/case03_firehose';
subcase = 'Astfalk2017Fig1_OFHI_pbk8'; 
dataFile = 'saveData_Astfalk2017Fig1_OFHI_pbk8.mat';
%
main_script = fullfile(case_path,subcase,'main_bopbk.m');
plot_script = fullfile(case_path,subcase,'selected_plot','bo_wpdat.m');
pbkdata = fullfile(case_path,subcase,'selected_plot', dataFile);
%
fprintf('\n========================================\n');
fprintf('Running: %s\n', main_script);
%
if exist(main_script, 'file') && exist(plot_script, 'file')
    run(main_script);
    run(plot_script);
else
    warning('Script not found: %s or %s', main_script, plot_script);
end

% Save results
currentPath = pwd;
datapath = './results/plt_Fig7/pbkData';
run('./tools/createDateFile(currentPath,datapath)');
copyfile(pbkdata,datapath);


%% Case#4.1 (currentDriven_IC: Tez/Tiz=2, bm)
clear; close all;
%
case_path = './examples/case04.1_currentDriven_IC';
subcase = 'IC_Fig6_TezTiz=2_bm'; 
dataFile = 'Basu2011PoPFig6_TezTiz=2_bm.mat';
%
main_script = fullfile(case_path,subcase,'main_bopbk.m');
plot_script = fullfile(case_path,subcase,'selected_plot','bo_wpdat.m');
pbkdata = fullfile(case_path,subcase,'selected_plot', dataFile);
%
fprintf('\n========================================\n');
fprintf('Running: %s\n', main_script);
%
if exist(main_script, 'file') && exist(plot_script, 'file')
    run(main_script);
    run(plot_script);
else
    warning('Script not found: %s or %s', main_script, plot_script);
end

% Save results
currentPath = pwd;
datapath = './results/plt_Fig9-10/pbkData_Fig9';
run('./tools/createDateFile(currentPath,datapath)');
copyfile(pbkdata,datapath);


%% Case#4.1 (currentDriven_IC: Tez/Tiz=2, km3)
clear; close all;
%
case_path = './examples/case04.1_currentDriven_IC';
subcase = 'IC_Fig6_TezTiz=2_km3'; 
dataFile = 'Basu2011PoPFig6_TezTiz=2_km3.mat';
%
main_script = fullfile(case_path,subcase,'main_bopbk.m');
plot_script = fullfile(case_path,subcase,'selected_plot','bo_wpdat.m');
pbkdata = fullfile(case_path,subcase,'selected_plot', dataFile);
%
fprintf('\n========================================\n');
fprintf('Running: %s\n', main_script);
%
if exist(main_script, 'file') && exist(plot_script, 'file')
    run(main_script);
    run(plot_script);
else
    warning('Script not found: %s or %s', main_script, plot_script);
end

% Save results
currentPath = pwd;
datapath = './results/plt_Fig9-10/pbkData_Fig9';
run('./tools/createDateFile(currentPath,datapath)');
copyfile(pbkdata,datapath);


%% Case#4.1 (currentDriven_IC: Tez/Tiz=5, bm)
clear; close all;
%
case_path = './examples/case04.1_currentDriven_IC';
subcase = 'IC_Fig6_TezTiz=5_bm'; 
dataFile = 'Basu2011PoPFig6_TezTiz=5_bm.mat';
%
main_script = fullfile(case_path,subcase,'main_bopbk.m');
plot_script = fullfile(case_path,subcase,'selected_plot','bo_wpdat.m');
pbkdata = fullfile(case_path,subcase,'selected_plot', dataFile);
%
fprintf('\n========================================\n');
fprintf('Running: %s\n', main_script);
%
if exist(main_script, 'file') && exist(plot_script, 'file')
    run(main_script);
    run(plot_script);
else
    warning('Script not found: %s or %s', main_script, plot_script);
end

% Save results
currentPath = pwd;
datapath = './results/plt_Fig9-10/pbkData_Fig9';
run('./tools/createDateFile(currentPath,datapath)');
copyfile(pbkdata,datapath);


%% Case#4.1 (currentDriven_IC: Tez/Tiz=5, km3)
clear; close all;
%
case_path = './examples/case04.1_currentDriven_IC';
subcase = 'IC_Fig6_TezTiz=5_km3'; 
dataFile = 'Basu2011PoPFig6_TezTiz=5_km3.mat';
%
main_script = fullfile(case_path,subcase,'main_bopbk.m');
plot_script = fullfile(case_path,subcase,'selected_plot','bo_wpdat.m');
pbkdata = fullfile(case_path,subcase,'selected_plot', dataFile);
%
fprintf('\n========================================\n');
fprintf('Running: %s\n', main_script);
%
if exist(main_script, 'file') && exist(plot_script, 'file')
    run(main_script);
    run(plot_script);
else
    warning('Script not found: %s or %s', main_script, plot_script);
end

% Save results
currentPath = pwd;
datapath = './results/plt_Fig9-10/pbkData_Fig9';
run('./tools/createDateFile(currentPath,datapath)');
copyfile(pbkdata,datapath);


%% Case#4.1 (currentDriven_IC: Tiz/Tix=2, bm)
clear; close all;
%
case_path = './examples/case04.1_currentDriven_IC';
subcase = 'IC_Fig7_TizTix=2_bm'; 
dataFile = 'Basu2011PoPFig7_TizTix=2_bm.mat';
%
main_script = fullfile(case_path,subcase,'main_bopbk.m');
plot_script = fullfile(case_path,subcase,'selected_plot','bo_wpdat.m');
pbkdata = fullfile(case_path,subcase,'selected_plot', dataFile);
%
fprintf('\n========================================\n');
fprintf('Running: %s\n', main_script);
%
if exist(main_script, 'file') && exist(plot_script, 'file')
    run(main_script);
    run(plot_script);
else
    warning('Script not found: %s or %s', main_script, plot_script);
end

% Save results
currentPath = pwd;
datapath = './results/plt_Fig9-10/pbkData_Fig10';
run('./tools/createDateFile(currentPath,datapath)');
copyfile(pbkdata,datapath);


%% Case#4.1 (currentDriven_IC: Tiz/Tix=2, km3)
clear; close all;
%
case_path = './examples/case04.1_currentDriven_IC';
subcase = 'IC_Fig7_TizTix=2_km3'; 
dataFile = 'Basu2011PoPFig7_TizTix=2_km3.mat';
%
main_script = fullfile(case_path,subcase,'main_bopbk.m');
plot_script = fullfile(case_path,subcase,'selected_plot','bo_wpdat.m');
pbkdata = fullfile(case_path,subcase,'selected_plot', dataFile);
%
fprintf('\n========================================\n');
fprintf('Running: %s\n', main_script);
%
if exist(main_script, 'file') && exist(plot_script, 'file')
    run(main_script);
    run(plot_script);
else
    warning('Script not found: %s or %s', main_script, plot_script);
end

% Save results
currentPath = pwd;
datapath = './results/plt_Fig9-10/pbkData_Fig10';
run('./tools/createDateFile(currentPath,datapath)');
copyfile(pbkdata,datapath);


%% Case#4.1 (currentDriven_IC: Tiz/Tix=5, bm)
clear; close all;
%
case_path = './examples/case04.1_currentDriven_IC';
subcase = 'IC_Fig7_TizTix=5_bm'; 
dataFile = 'Basu2011PoPFig7_TizTix=5_bm.mat';
%
main_script = fullfile(case_path,subcase,'main_bopbk.m');
plot_script = fullfile(case_path,subcase,'selected_plot','bo_wpdat.m');
pbkdata = fullfile(case_path,subcase,'selected_plot', dataFile);
%
fprintf('\n========================================\n');
fprintf('Running: %s\n', main_script);
%
if exist(main_script, 'file') && exist(plot_script, 'file')
    run(main_script);
    run(plot_script);
else
    warning('Script not found: %s or %s', main_script, plot_script);
end

% Save results
currentPath = pwd;
datapath = './results/plt_Fig9-10/pbkData_Fig10';
run('./tools/createDateFile(currentPath,datapath)');
copyfile(pbkdata,datapath);

% Save results for Fig11
datapath = './results/plt_Fig11/pbkData';
run('./tools/createDateFile(currentPath,datapath)');
copyfile(pbkdata,datapath);


%% Case#4.1 (currentDriven_IC: Tiz/Tix=5, km3)
clear; close all;
%
case_path = './examples/case04.1_currentDriven_IC';
subcase = 'IC_Fig7_TizTix=5_km3'; 
dataFile = 'Basu2011PoPFig7_TizTix=5_km3.mat';
%
main_script = fullfile(case_path,subcase,'main_bopbk.m');
plot_script = fullfile(case_path,subcase,'selected_plot','bo_wpdat.m');
pbkdata = fullfile(case_path,subcase,'selected_plot', dataFile);
%
fprintf('\n========================================\n');
fprintf('Running: %s\n', main_script);
%
if exist(main_script, 'file') && exist(plot_script, 'file')
    run(main_script);
    run(plot_script);
else
    warning('Script not found: %s or %s', main_script, plot_script);
end

% Save results
currentPath = pwd;
datapath = './results/plt_Fig9-10/pbkData_Fig10';
run('./tools/createDateFile(currentPath,datapath)');
copyfile(pbkdata,datapath);
%
% Save results for Fig11
datapath = './results/plt_Fig11/pbkData';
run('./tools/createDateFile(currentPath,datapath)');
copyfile(pbkdata,datapath);

%% Case#4.1 (currentDriven_IC: Tiz/Tix=5, sigma=0.5, bm)
clear; close all;
%
case_path = './examples/case04.1_currentDriven_IC';
subcase = 'IC_Fig7_TizTix=5_sig=0.5_bm'; 
dataFile = 'Basu2011PoPFig7_TizTix=5_sig=0.5_bm.mat';
%
main_script = fullfile(case_path,subcase,'main_bopbk.m');
plot_script = fullfile(case_path,subcase,'selected_plot','bo_wpdat.m');
pbkdata = fullfile(case_path,subcase,'selected_plot', dataFile);
%
fprintf('\n========================================\n');
fprintf('Running: %s\n', main_script);
%
if exist(main_script, 'file') && exist(plot_script, 'file')
    run(main_script);
    run(plot_script);
else
    warning('Script not found: %s or %s', main_script, plot_script);
end

% Save results
currentPath = pwd;
datapath = './results/plt_Fig11/pbkData';
run('./tools/createDateFile(currentPath,datapath)');
copyfile(pbkdata,datapath);

%% Case#4.1 (currentDriven_IC: Tiz/Tix=5, sigma=0.5, km3)
clear; close all;
%
case_path = './examples/case04.1_currentDriven_IC';
subcase = 'IC_Fig7_TizTix=5_sig=0.5_km3'; 
dataFile = 'Basu2011PoPFig7_TizTix=5_sig=0.5_km3.mat';
%
main_script = fullfile(case_path,subcase,'main_bopbk.m');
plot_script = fullfile(case_path,subcase,'selected_plot','bo_wpdat.m');
pbkdata = fullfile(case_path,subcase,'selected_plot', dataFile);
%
fprintf('\n========================================\n');
fprintf('Running: %s\n', main_script);
%
if exist(main_script, 'file') && exist(plot_script, 'file')
    run(main_script);
    run(plot_script);
else
    warning('Script not found: %s or %s', main_script, plot_script);
end

% Save results
currentPath = pwd;
datapath = './results/plt_Fig11/pbkData';
run('./tools/createDateFile(currentPath,datapath)');
copyfile(pbkdata,datapath);


%% Case#4.1 (currentDriven_IC: Tiz/Tix=5, sigma=1, bm)
clear; close all;
%
case_path = './examples/case04.1_currentDriven_IC';
subcase = 'IC_Fig7_TizTix=5_sig=1_bm'; 
dataFile = 'Basu2011PoPFig7_TizTix=5_sig=1_bm.mat';
%
main_script = fullfile(case_path,subcase,'main_bopbk.m');
plot_script = fullfile(case_path,subcase,'selected_plot','bo_wpdat.m');
pbkdata = fullfile(case_path,subcase,'selected_plot', dataFile);
%
fprintf('\n========================================\n');
fprintf('Running: %s\n', main_script);
%
if exist(main_script, 'file') && exist(plot_script, 'file')
    run(main_script);
    run(plot_script);
else
    warning('Script not found: %s or %s', main_script, plot_script);
end

% Save results
currentPath = pwd;
datapath = './results/plt_Fig11/pbkData';
run('./tools/createDateFile(currentPath,datapath)');
copyfile(pbkdata,datapath);


%% Case#4.1 (currentDriven_IC: Tiz/Tix=5, sigma=1, km3)
clear; close all;
%
case_path = './examples/case04.1_currentDriven_IC';
subcase = 'IC_Fig7_TizTix=5_sig=1_km3'; 
dataFile = 'Basu2011PoPFig7_TizTix=5_sig=1_km3.mat';
%
main_script = fullfile(case_path,subcase,'main_bopbk.m');
plot_script = fullfile(case_path,subcase,'selected_plot','bo_wpdat.m');
pbkdata = fullfile(case_path,subcase,'selected_plot', dataFile);
%
fprintf('\n========================================\n');
fprintf('Running: %s\n', main_script);
%
if exist(main_script, 'file') && exist(plot_script, 'file')
    run(main_script);
    run(plot_script);
else
    warning('Script not found: %s or %s', main_script, plot_script);
end

% Save results
currentPath = pwd;
datapath = './results/plt_Fig11/pbkData';
run('./tools/createDateFile(currentPath,datapath)');
copyfile(pbkdata,datapath);


%% Case#4.2 (EMIC, theta=0)
clear; close all;
%
case_path = './examples/case04.2_EMIC';
subcase = 'bai2025CPCfig13_EMIC_th=0_km1'; 
dataFile = 'pbk_EMIC_th=0_B0=100nT_Tsz=1879.mat';
%
main_script = fullfile(case_path,subcase,'main_bopbk.m');
plot_script = fullfile(case_path,subcase,'selected_plot','bo_wpdat.m');
pbkdata = fullfile(case_path,subcase,'selected_plot', dataFile);
%
fprintf('\n========================================\n');
fprintf('Running: %s\n', main_script);
%
if exist(main_script, 'file') && exist(plot_script, 'file')
    run(main_script);
    run(plot_script);
else
    warning('Script not found: %s or %s', main_script, plot_script);
end

% Save results
currentPath = pwd;
datapath = './results/plt_Fig12/pbkData';
run('./tools/createDateFile(currentPath,datapath)');
copyfile(pbkdata,datapath);


%% Case#4.2 (EMIC, theta=15)
clear; close all;
%
case_path = './examples/case04.2_EMIC';
subcase = 'bai2025CPCfig13_EMIC_th=15_km1'; 
dataFile = 'pbk_EMIC_th=15_B0=100nT_Tsz=1879.mat';
%
main_script = fullfile(case_path,subcase,'main_bopbk.m');
plot_script = fullfile(case_path,subcase,'selected_plot','bo_wpdat.m');
pbkdata = fullfile(case_path,subcase,'selected_plot', dataFile);
%
fprintf('\n========================================\n');
fprintf('Running: %s\n', main_script);
%
if exist(main_script, 'file') && exist(plot_script, 'file')
    run(main_script);
    run(plot_script);
else
    warning('Script not found: %s or %s', main_script, plot_script);
end

% Save results
currentPath = pwd;
datapath = './results/plt_Fig12/pbkData';
run('./tools/createDateFile(currentPath,datapath)');
copyfile(pbkdata,datapath);



%% Case#4.2 (EMIC, theta=40)
clear; close all;
%
case_path = './examples/case04.2_EMIC';
subcase = 'bai2025CPCfig13_EMIC_th=40_km1'; 
dataFile = 'pbk_EMIC_th=40_B0=100nT_Tsz=1879.mat';
%
main_script = fullfile(case_path,subcase,'main_bopbk.m');
plot_script = fullfile(case_path,subcase,'selected_plot','bo_wpdat.m');
pbkdata = fullfile(case_path,subcase,'selected_plot', dataFile);
%
fprintf('\n========================================\n');
fprintf('Running: %s\n', main_script);
%
if exist(main_script, 'file') && exist(plot_script, 'file')
    run(main_script);
    run(plot_script);
else
    warning('Script not found: %s or %s', main_script, plot_script);
end

% Save results
currentPath = pwd;
datapath = './results/plt_Fig12/pbkData';
run('./tools/createDateFile(currentPath,datapath)');
copyfile(pbkdata,datapath);

%%
% Case#4.2 (EMIC, theta=55)
clear; close all;
%
case_path = './examples/case04.2_EMIC';
subcase = 'bai2025CPCfig13_EMIC_th=55_km1'; 
dataFile = 'pbk_EMIC_th=55_B0=100nT_Tsz=1879.mat';
%
main_script = fullfile(case_path,subcase,'main_bopbk.m');
plot_script = fullfile(case_path,subcase,'selected_plot','bo_wpdat.m');
pbkdata = fullfile(case_path,subcase,'selected_plot', dataFile);
%
fprintf('\n========================================\n');
fprintf('Running: %s\n', main_script);
%
if exist(main_script, 'file') && exist(plot_script, 'file')
    run(main_script);
    run(plot_script);
else
    warning('Script not found: %s or %s', main_script, plot_script);
end

% Save results
currentPath = pwd;
datapath = './results/plt_Fig12/pbkData';
run('./tools/createDateFile(currentPath,datapath)');
copyfile(pbkdata,datapath);


%% Case#3 (OFHI_bm 3D; OFHI_pbk 3D)
clear; close all;
%
% run_case3_firehose3D;


%% Case#4.2 (EMIC 3D)
clear; close all;
%
% run_case4dot2_EMIC3D;

clear; close all;


% ================== LOG COMPLETION ====================
fprintf('Done at: %s\n', datetime('now'));

% Stop logging
diary off;
rmpath('./tools');
