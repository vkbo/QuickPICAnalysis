
%
%  Class Object :: Wrapper for QuickPIC data sets
% ************************************************
%

classdef QPICData
    
    properties(GetAccess='public', SetAccess='public')

        Path        = '';    % Path to dataset
        PathID      = '';    % Path as ID instead of free text input
        Config      = [];    % Content of the config files and extraction of all runtime variables
        Silent      = false; % Set to 1 to disable command window output

    end % properties

    properties(GetAccess='public', SetAccess='private')

        Elements    = {};    % Struct of all datafiles in dataset
        MSData      = {};    % Struct of all MS data
        DataSets    = {};    % Available datasets in folders indicated by LocalConfig.m
        DefaultPath = {};    % Default data folder
        Temp        = '';    % Temp folder (set in LocalConfig.m)
        HasData     = false; % True if any data folder exists
        Consistent  = false; % True if all data folders have the same number of files

    end % properties

    properties(GetAccess='private', SetAccess='private')

        DefaultData = {};    % Data in default folder

    end % properties

    %
    % Constructor
    %
    
    methods
        
        function obj = QPICData(varargin)
            
            % Parse input
            oOpt = inputParser;
            addParameter(oOpt, 'Silent', 'No');
            parse(oOpt, varargin{:});
            stOpt = oOpt.Results;
            
            if strcmpi(stOpt.Silent, 'Yes')
                obj.Silent = true;
            end % if

            % Initiate QPICData
            LocalConfig;
            
            obj.Temp          = sLocalTemp;
            obj.Config        = QPICConfig;
            obj.Config.Silent = obj.Silent;
            
            obj.DefaultPath = stFolders;
            if ~obj.Silent
                fprintf('Scanning default data folder(s)\n');
            end % if
            
            stFields = fieldnames(obj.DefaultPath);

            for f=1:length(stFields)

                sName  = stFields{f};
                sPath  = obj.DefaultPath.(stFields{f}).Path;
                iDepth = obj.DefaultPath.(stFields{f}).Depth;
                
                obj.DefaultPath.(stFields{f}).Available = 0;

                if isdir(sPath)
                    
                    if ~obj.Silent
                        fprintf('Scanning %s\n', sPath);
                    end % if
                    obj.DefaultPath.(stFields{f}).Available = 1;
                    
                    stScan.(sName)(1) = struct('Path', sPath, 'Name', sName, 'Level', 0);
                    for r=0:iDepth-1
                        for p=1:length(stScan.(sName))
                            if stScan.(sName)(p).Level == r
                                stDir = dir(stScan.(sName)(p).Path);
                                for d=1:length(stDir)
                                    if stDir(d).isdir == 1 && ~strcmp(stDir(d).name, '.') && ~strcmp(stDir(d).name, '..')
                                        stScan.(sName)(end+1) = struct('Path',  [stScan.(sName)(p).Path, '/', stDir(d).name], ...
                                                                       'Name',  stDir(d).name, ...
                                                                       'Level', r+1);
                                    end % if
                                end % for
                            end % if
                        end % for
                    end % for
                    
                    for s=1:length(stScan.(sName))
                        if stScan.(sName)(s).Level == iDepth
                            sSet = structname(stScan.(sName)(s).Name);

                            obj.DataSets.ByName.(sSet).Path      = stScan.(sName)(s).Path;
                            obj.DataSets.ByName.(sSet).Name      = stScan.(sName)(s).Name;
                            obj.DataSets.ByName.(sSet).HasData   = isdir([stScan.(sName)(s).Path, '/MS']);

                            obj.DataSets.ByPath.(sName).(sSet).Path      = stScan.(sName)(s).Path;
                            obj.DataSets.ByPath.(sName).(sSet).Name      = stScan.(sName)(s).Name;
                            obj.DataSets.ByPath.(sName).(sSet).HasData   = isdir([stScan.(sName)(s).Path, '/MS']);
                        end % if
                    end % for
                end % if
            end % for

        end % function

    end % methods
    
end % classdef
