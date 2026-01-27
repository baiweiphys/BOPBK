
clear; close all;

%% Case#4.2 (EMIC 3D)
case_path = './examples/case04.2_EMIC';
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
datapath = './results/plt_Fig13/pbkData';
run('./tools/createDateFile(currentPath,datapath)');
copyfile(pbkdata,datapath);

clear; close all;