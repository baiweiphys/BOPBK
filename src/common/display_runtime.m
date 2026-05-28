function display_runtime(runtime)
% DISPLAY_RUNTIME  Prints the execution time.
% Input: runtime - Total elapsed time in seconds (usually from toc).
% Author: Bai Wei (baiweiphys@gmail.com)
% Date: 2026.01.30

    if runtime < 60
        % Case 1: Less than 1 minute - display seconds with 2 decimal places
        fprintf('Elapsed time: %.2f s\n', runtime);
        
    elseif runtime < 3600
        % Case 2: Less than 1 hour - display minutes and seconds
        m = floor(runtime / 60);
        s = mod(runtime, 60);
        fprintf('Elapsed time: %d m %d s\n', m, round(s));
        
    else
        % Case 3: 1 hour or more - display hours, minutes, and seconds
        h = floor(runtime / 3600);
        m = floor(mod(runtime, 3600) / 60);
        s = mod(runtime, 60);
        fprintf('Elapsed time: %d h %d m %d s\n', h, m, round(s));
    end
end