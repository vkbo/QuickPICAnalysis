
%
%  Local Settings
% ****************
%  Copy this file as LocalConfig.m and change the settings below
%

% Local temp directory
sLocalTemp = '/home/vkbo/Temp';

%  List of data folders to scan
% *******************************
%  Depth indicates at which sub level the data folders are stored.
%  If the same dataset exists in several locations, the last entry in the
%  list takes priority.

[~,sHost] = system('hostname');
sHost     = strtrim(sHost);

% Rey
if strcmpi(sHost,'Rey')
    stFolders.IntTSeries = struct('Path', '/scratch/QuickPICData/T-Series', 'Depth', 1, 'Name', 'Data: T-Series');
    stFolders.IntTSeries = struct('Path', '/scratch/SimData/Q-Series',      'Depth', 1, 'Name', 'Data: Q-Series');
end % if
