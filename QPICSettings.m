
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

% Trinity
if strcmpi(sHost,'Trinity')
    stFolders.IntTSeries = struct('Path', '/scratch/QuickPICData/T-Series', 'Depth', 1, 'Name', 'Data: T-Series');
    stFolders.IntQSeries = struct('Path', '/data/SimData/Q-Series',         'Depth', 1, 'Name', 'Data: Q-Series');
end % if

% Leela
if strcmpi(sHost,'Leela')
    stFolders.IntQSeries = struct('Path', '/scratch/SimData/Q-Series',         'Depth', 1, 'Name', 'Data: Q-Series');
end % if

% Xena
if strcmpi(sHost,'Xena')
    stFolders.ExtQSeries = struct('Path', '/media/vkbo/DataDrive/SimData/Q-Series', 'Depth', 1, 'Name', 'Ext: Q-Series');
    stFolders.IntQSeries = struct('Path', '/scratch/SimData/Q-Series',              'Depth', 1, 'Name', 'Data: Q-Series');
end % if
