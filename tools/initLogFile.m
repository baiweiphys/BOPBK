function logFilePath = initLogFile(scriptName)
% INITLOGFILE Create timestamped log file
%   INPUT:  scriptName - Name of the calling script (e.g., 'plotall')
%   OUTPUT: logFilePath - Full path to the created log file
    
    % Create logs directory if it doesn't exist
    logDir = 'logs';
    if ~exist(logDir, 'dir')
        mkdir(logDir);
        fprintf('Created log directory: %s\n', logDir);
    end
    
    % Generate filename with timestamp
    timeStr = datestr(now, 'yyyymmdd_HHMMSS');  % Format: YYYYMMDD_HHMMSS
    filename = sprintf('%s_%s.txt', scriptName, timeStr);
    
    % Construct full file path
    logFilePath = fullfile(logDir, filename);
    
    % Ensure file is new (optional - removes existing file)
    if exist(logFilePath, 'file')
        delete(logFilePath);
    end
    
    fprintf('Log initialized: %s\n', logFilePath);
end
