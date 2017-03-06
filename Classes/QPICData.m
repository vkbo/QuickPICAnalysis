
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
        SimData     = {};    % Struct of all simulation data
        DataSets    = {};    % Available datasets in folders indicated by QPICSettings.m
        DefaultPath = {};    % Default data folder
        Temp        = '';    % Temp folder (set in QPICSettings.m)
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
            QPICSettings;
            
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
                            
                            allFiles = dir(stScan.(sName)(s).Path);
                            allDirs  = allFiles([allFiles(:).isdir]);
                            numDirs  = numel(allDirs);
                            iHasData = numDirs > 3;

                            obj.DataSets.ByName.(sSet).Path            = stScan.(sName)(s).Path;
                            obj.DataSets.ByName.(sSet).Name            = stScan.(sName)(s).Name;
                            obj.DataSets.ByName.(sSet).HasData         = iHasData;

                            obj.DataSets.ByPath.(sName).(sSet).Path    = stScan.(sName)(s).Path;
                            obj.DataSets.ByPath.(sName).(sSet).Name    = stScan.(sName)(s).Name;
                            obj.DataSets.ByPath.(sName).(sSet).HasData = iHasData;
                        end % if
                    end % for
                end % if
            end % for

        end % function

    end % methods
    
    %
    % Setters and Getters
    %

    methods
        
        function obj = set.Path(obj, vInput)
            
            %
            %  Sets obj.Path and scans data tree
            % ***********************************
            %
            
            iHasData = 0;
            
            if isstruct(vInput)
                if isfield(vInput, 'Path')
                    obj.Path = vInput.Path;
                    iHasData = vInput.HasData;
                else
                    fprintf(2, 'Error: Path not recognised or found.\n');
                    return;
                end % if
            elseif isdir(vInput)
                obj.Path = vInput;
                allFiles = dir(stScan.(sName)(s).Path);
                allDirs  = allFiles([allFiles(:).isdir]);
                numDirs  = numel(allDirs);
                iHasData = numDirs > 3;
            else
                sField = structname(vInput);
                if isfield(obj.DataSets.ByName, sField)
                    obj.Path = obj.DataSets.ByName.(sField).Path;
                    iHasData = obj.DataSets.ByName.(sField).HasData;
                else
                    fprintf(2, 'Error: Path not recognised or found.\n');
                    return;
                end % if
            end % if

            if ~obj.Silent
                fprintf('Path is %s\n', obj.Path);
            end % if
            
            % Scanning simulation folder
            obj.Elements = obj.fScanFolder;
            obj.SimData  = obj.fScanElements;

            % Set path in QPICConfig object
            obj.Config.Path = obj.Path;
            obj.HasData     = iHasData;
            
            if mod(obj.SimData.MaxFiles,obj.SimData.MinFiles) == 0
               obj.Consistent = true;
            else
               obj.Consistent = false;
            end % if
            
            %obj.Translate = Variables(obj.Config.Simulation.Coordinates, obj.RunningZ);
            
            % Output Dataset Info
            if ~obj.Silent
                if obj.HasData
                    fprintf('Folder contains simulation data.\n');
                end % if
                if ~obj.Consistent
                    fprintf('Simulation has varying number of time dumps.\n');
                end % if
            end % if

        end % function
        
        function obj = set.Elements(obj, stElements)

            obj.Elements = stElements;

        end % function
        
    end % methods

    %
    % Public Methods
    %
    
    methods(Access = 'public')
        
        function Version(~)
            
            fprintf('QPICAnalysis Version Dev0.1\n');
            
        end % function

        function aReturn = Data(obj, iTime, sType, sSet, sSpecies, sSlice)
            
            %
            %  Data-extraction function
            % **************************
            %
            %  Input:
            % ========
            %  iTime    :: Time dump to extract
            %  sType    :: Data type [F, J, PSI, Q]
            %  sSet     :: Data set i.e. EX, BY, etc
            %  sSpecies :: Particle species i.e. EB, EP01, etc
            %  sSlice   :: Which data slice i.e. XY, XZ or YZ or blank for 3D
            %
            
            % Input/Output
            aReturn = [];

            if nargin == 1
                fprintf('\n');
                fprintf('  Data-extraction function\n');
                fprintf(' **************************\n');
                fprintf('\n');
                fprintf('  Input:\n');
                fprintf(' ========\n');
                fprintf('  iTime    :: Time dump to extract\n');
                fprintf('  sType    :: Data type [F, J, PSI, Q]\n');
                fprintf('  sSet     :: Data set i.e. EX, BY, etc\n');
                fprintf('  sSpecies :: Particle species i.e. EB, EP01, etc\n');
                fprintf('  sSlice   :: Which data slice i.e. XY, XZ or YZ or blank for 3D\n');
                fprintf('\n');
                return;
            end % if
            
            if strcmp(obj.Path, '')
                fprintf(2, 'Error: No dataset has been loaded.\n');
                return;
            end % if
            
            % Enforce upper case
            sType    = upper(sType);
            sSet     = upper(sSet);
            sSpecies = upper(sSpecies);
            sSlice   = upper(sSlice);
            
            if isempty(sSlice)
                sSlice = 'All';
            end % if
            if isempty(sSet)
                sSet = 'None';
            end % if
            if isempty(sSpecies)
                sSpecies = 'None';
            end % if

            if isempty(sType)
                fprintf(2, 'Error: Data type needs to be specified.\n');
                return;
            end % if
            if ~obj.DataSetExists(sType, sSet, sSpecies, sSlice)
                fprintf(2, 'Error: Specified data set does not exist.\n');
                return;
            end % if

            % Extract path
            iIndex    = obj.SimData.Index.(sType).(sSet).(sSpecies).(sSlice);
            sFolder   = obj.SimData.Data(iIndex).Path;
            iFiles    = obj.SimData.Data(iIndex).Files;
            sTimeNExt = sprintf('%08d', iTime);
            sFile     = ['/',sFolder,'_',sTimeNExt,'.h5'];
            sLoad     = [obj.Path,'/',sFolder,sFile];
            
            % Check if datafile exists
            if iTime > iFiles
                fprintf(2, 'Error: Dump %d does not exist. Last dump is %d.\n', iTime, iFiles);
                return;
            end % if

            % Wraps reading in a try/catch because sometimes files are corrupt even if they exist.
            try
                aReturn = double(h5read(sLoad, ['/',sFolder]));
            catch
                fprintf(2, 'Error reading file %s\n', sLoad);
            end % try
            
        end % function
        
        function bReturn = DataSetExists(obj, sType, sSet, sSpecies, sSlice)
            
            bReturn = false;
            
            if isempty(sSlice)
                sSlice = 'All';
            end % if

            [~,iMS] = size(obj.SimData.Data);
            for m=1:iMS
                if   strcmpi(obj.SimData.Data(m).Type, sType) ...
                  && strcmpi(obj.SimData.Data(m).Set, sSet) ...
                  && strcmpi(obj.SimData.Data(m).Species, sSpecies) ...
                  && strcmpi(obj.SimData.Data(m).Slice, sSlice)
                    bReturn = true;
                    return;
                end % if
            end % for
            
        end % function

    end % methods
    
    %
    % Private Methods
    %
    
    methods(Access = 'private')
        
        function stReturn = fScanFolder(obj)
            
            stReturn.Info = struct('Path', '', 'Dirs', 0, 'Files', 0);

            allFiles = dir(obj.Path);
            allDirs  = allFiles([allFiles(:).isdir]);
            cExclude = {'.','..','ELOG'};
            
            for f=1:numel(allDirs)
                
                sPath = allDirs(f).name;

                if sum(ismember(sPath, cExclude)) > 0
                    continue;
                end % if
                
                cElems = strsplit(sPath, '-');
                sName  = cElems{1};
                if numel(cElems) > 1
                    sDir = cElems{2};
                else
                    sDir = 'All'; % If there are no directions specified, assume it's for all directions
                end % if

                subFiles = dir([obj.Path '/' sPath]);
                numFiles = numel(subFiles) - 2;

                stReturn.(sName).Info        = struct('Path', '', 'Dirs', 0, 'Files', 0);
                stReturn.(sName).(sDir).Info = struct('Path', sPath, 'Dirs', 0, 'Files', numFiles);
                
            end % for
            
            cDataSets = fieldnames(stReturn);
            stReturn.Info.Dirs = numel(cDataSets) - 1;
            
            for s=1:numel(cDataSets)

                sName = cDataSets{s};
                
                if strcmpi(sName,'Info')
                    continue;
                end % if
                
                cDirections = fieldnames(stReturn.(sName));
                stReturn.(sName).Info.Dirs = numel(cDirections) - 1;

            end % for
        
        end % function
        
        function stReturn = fScanElements(obj)
            
            stReturn = struct();
            stData   = struct();
            stIndex  = struct();
            iRow     = 1;
            iMax     = 0;
            iMin     = 1e6;
            
            stType = fieldnames(obj.Elements);
            for i=2:length(stType)
                
                sType     = stType{i};
                sDType    = '';
                sDSet     = '';
                sDSpecies = '';
                
                switch(sType(1))
                    case 'F'
                        sDType    = 'F';
                        sDSet     = sType(2:end);
                        sDSpecies = 'None';
                    case 'J'
                        sDType    = 'J';
                        sDSet     = ['J' sType(end)];
                        sDSpecies = 'None';
                    case 'P'
                        sDType    = 'PSI';
                        sDSet     = 'None';
                        sDSpecies = 'None';
                    case 'Q'
                        sDType    = 'Q';
                        sDSet     = 'None';
                        sDSpecies = sType(2:end);
                end % switch
                
                if ~isempty(sDType)
                    stSlice = fieldnames(obj.Elements.(sType));
                    for j=2:length(stSlice)
                        stData(iRow).Type    = sDType;
                        stData(iRow).Set     = sDSet;
                        stData(iRow).Species = sDSpecies;
                        stData(iRow).Slice   = stSlice{j};
                        stData(iRow).Path    = obj.Elements.(sType).(stSlice{j}).Info.Path;
                        stData(iRow).Files   = obj.Elements.(sType).(stSlice{j}).Info.Files;
                        stIndex.(sDType).(sDSet).(sDSpecies).(stSlice{j}) = iRow;
                        
                        iFiles = stData(iRow).Files;
                        if iFiles > iMax
                            iMax = iFiles;
                        end % if
                        if iFiles < iMin
                            iMin = iFiles;
                        end % if
                        iRow = iRow + 1;
                    end % for
                end % if

            end % for
            
            if iMin == 1e6
                iMin = 0;
            end % if
            
            stReturn.Data     = stData;
            stReturn.Index    = stIndex;
            stReturn.MinFiles = iMin;
            stReturn.MaxFiles = iMax;
            
        end % function

    end % methods

end % classdef
