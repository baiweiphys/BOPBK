
clear; close all;

%% Case#3 (OFHI_bm 3D)
case_path = './examples/case03_firehose';
subcase = 'Astfalk2017Fig1_OFHI3d_bm'; 
dataFile = 'output/Astfalk2017Fig1_OFHI3d_bm.mat';
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
datapath = './results/plt_Fig8/pbkData';
run('./tools/createDateFile(currentPath,datapath)');
copyfile(pbkdata,datapath);

%% Case#3 (OFHI_pbk210 3D)
clear; close all;
%
case_path = './examples/case03_firehose';
subcase = 'Astfalk2017Fig1_OFHI3d_pbk210'; 
dataFile = 'output/Astfalk2017Fig1_OFHI3d_pbk210.mat';
%
main_script = fullfile(case_path,subcase,'main_bopbk.m');
pbkdata = fullfile(case_path,subcase,dataFile);
%
fprintf('\n========================================\n');
fprintf('Running: %s\n', main_script);
fprintf('========================================\n');
%
if exist(main_script, 'file')
    run(main_script);
else
    warning('Script not found: %s', main_script);
end

% Save results
currentPath = pwd;
datapath = './results/plt_Fig8/pbkData';
run('./tools/createDateFile(currentPath,datapath)');
copyfile(pbkdata,datapath);

%%
clear; close all;
