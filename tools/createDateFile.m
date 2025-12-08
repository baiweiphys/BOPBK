function createDateFile(currentPath,fname)
% @Description: Create a new data folder if it doesn't exist.
% @Filename: createDateFile.m
% @Author: Bai Wei (baiweiphys@gmail.com)
% @Date: 2023-10-01
% @LastEditors: Bai Wei
% @LastEditTime: 2023-11-15

datafolder = [currentPath,'/',fname];
if exist(datafolder)==0
    mkdir(datafolder);
else
    disp('The data folder exists.');
end