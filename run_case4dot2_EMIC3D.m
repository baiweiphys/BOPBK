% RUN_CASE4DOT2_EMIC3D Execute the EMIC 3D benchmark case.
% This script runs the configured case under examples/case04.2_EMIC
% and copies the generated output file into the corresponding results folder.

clear; close all;

%% Case#4.2 (EMIC 3D)
case_path = './examples/case04.3_EMIC_bai2025CPC';
subcase = 'bai2025CPCfig12_EMIC_allAngles'; 
dataFile = 'output/bopbkData.mat';
%
main_script = fullfile(case_path,subcase,'main_bopbk.m');
pbkdata = fullfile(case_path,subcase,dataFile);
%
fprintf('\n========================================\n');
fprintf('Running: %s\n', main_script);
%
if exist(main_script, 'file')
    run(main_script);
else
    warning('Script not found: %s', main_script);
end

% Save results
currentPath = pwd;
datapath = './results/plt_Fig15/pbkData';
run('./tools/createDateFile(currentPath,datapath)');
copyfile(pbkdata,datapath);

clear; close all;
