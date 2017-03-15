
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
            
            fprintf('QPICAnalysis Version Dev0.2\n');
            
        end % function

        function stReturn = BeamInfo(obj, iBeam)
            
            %
            %  Attempts to Calculate Beam Info
            % *********************************
            %
            
            stReturn = {};
            
            dC    = obj.Config.Constants.SI.SpeedOfLight;
            dE    = obj.Config.Constants.SI.ElementaryCharge;
            dLFac = obj.Config.Convert.SI.LengthFac;
            sBeam = sprintf('Beam%02d',iBeam);
            
            if ~obj.Silent
                fprintf('\n');
                fprintf(' Beam Info for %s\n',sBeam);
                fprintf('**********************\n');
                fprintf('\n');
            end % if

            %
            % Sigma and Mean
            %

            dSIMeanX1  = obj.Config.Beam.(sBeam).Position(1);
            dSIMeanX2  = obj.Config.Beam.(sBeam).Position(2);
            dSIMeanX3  = obj.Config.Beam.(sBeam).Position(3);
            
            dSISigmaX1 = obj.Config.Beam.(sBeam).Profile(1);
            dSISigmaX2 = obj.Config.Beam.(sBeam).Profile(2);
            dSISigmaX3 = obj.Config.Beam.(sBeam).Profile(3);

            dSIEmittX1 = obj.Config.Beam.(sBeam).Emittance(1);
            dSIEmittX2 = obj.Config.Beam.(sBeam).Emittance(2);
            dSIEmittX3 = obj.Config.Beam.(sBeam).Emittance(3);

            dBeamGamma   = obj.Config.Beam.(sBeam).Gamma;
            dBeamBeta    = sqrt(1-1/dBeamGamma^2);
            dBeamCharge  = obj.Config.Beam.(sBeam).BeamCharge;
            dBeamNum     = obj.Config.Beam.(sBeam).NumParticles;
            dBeamCurrent = dBeamCharge*dC*dBeamBeta / sqrt(2*pi*dSISigmaX1^2);

            stReturn.X1Mean      = dSIMeanX1;
            stReturn.X2Mean      = dSIMeanX2;
            stReturn.X3Mean      = dSIMeanX3;

            stReturn.X1Sigma     = dSISigmaX1;
            stReturn.X2Sigma     = dSISigmaX2;
            stReturn.X3Sigma     = dSISigmaX3;

            stReturn.X1Emittance = dSIEmittX1;
            stReturn.X2Emittance = dSIEmittX2;
            stReturn.X3Emittance = dSIEmittX3;

            stReturn.Particles   = dBeamNum;
            stReturn.Charge      = dBeamCharge;
            stReturn.Current     = dBeamCurrent;

            [dSIMeanX1,  sUnitM1] = fAutoScale(dSIMeanX1,'m',1e-9);
            [dSIMeanX2,  sUnitM2] = fAutoScale(dSIMeanX2,'m',1e-9);
            [dSIMeanX3,  sUnitM3] = fAutoScale(dSIMeanX3,'m',1e-9);

            [dSISigmaX1, sUnitS1] = fAutoScale(dSISigmaX1,'m',1e-9);
            [dSISigmaX2, sUnitS2] = fAutoScale(dSISigmaX2,'m',1e-9);
            [dSISigmaX3, sUnitS3] = fAutoScale(dSISigmaX3,'m',1e-9);

            [dSIEmittX1, sUnitE1] = fAutoScale(dSIEmittX1,'m',1e-9);
            [dSIEmittX2, sUnitE2] = fAutoScale(dSIEmittX2,'m',1e-9);
            [dSIEmittX3, sUnitE3] = fAutoScale(dSIEmittX3,'m',1e-9);
            
            [dBeamGamma,   sGammaUnit]   = fAutoScale(dBeamGamma,'');
            [dBeamCharge,  sChargeUnit]  = fAutoScale(dBeamCharge,'C');
            [dBeamCurrent, sCurrentUnit] = fAutoScale(dBeamCurrent,'A');
            [dBeamNum,     sBeamNumUnit] = fAutoScale(dBeamNum,'');

            if ~obj.Silent
                fprintf(' X1 Mean, Sigma: %8.3f %2s, %8.3f %2s\n',dSIMeanX1,sUnitM1,dSISigmaX1,sUnitS1);
                fprintf(' X2 Mean, Sigma: %8.3f %2s, %8.3f %2s\n',dSIMeanX2,sUnitM2,dSISigmaX2,sUnitS2);
                fprintf(' X3 Mean, Sigma: %8.3f %2s, %8.3f %2s\n',dSIMeanX3,sUnitM3,dSISigmaX3,sUnitS3);
                fprintf('\n');
                fprintf(' X1 Emittance:   %8.3f %2s\n',dSIEmittX1,sUnitE1);
                fprintf(' X2 Emittance:   %8.3f %2s\n',dSIEmittX2,sUnitE2);
                fprintf(' X3 Emittance:   %8.3f %2s\n',dSIEmittX3,sUnitE3);
                fprintf('\n');
                fprintf(' Beam Gamma:     %8.3f %s\n', dBeamGamma, sGammaUnit);
                fprintf(' Beam Charge:    %8.3f %s\n', dBeamCharge, sChargeUnit);
                fprintf(' Beam Current:   %8.3f %s\n', dBeamCurrent, sCurrentUnit);
                fprintf(' Particle Count: %8.3f %s Particles\n', dBeamNum, sBeamNumUnit);
                fprintf('\n');
            end % if
            
        end % function

        function aReturn = Data(obj, iTime, sType, sSet, sSpecies, sSlice)
            
            %
            %  Data-extraction function
            % **************************
            %
            %  Input:
            % ========
            %  iTime    :: Time dump to extract
            %  sType    :: Data type i.e F, J, PSI, Q, RAW
            %  sSet     :: Data set i.e. EX, BY, PSI, etc
            %  sSpecies :: Particle species i.e. EB, EP01, etc or for RAW, 01, 02, etc
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
                fprintf('  sType    :: Data type i.e. F, J, PSI, Q, RAW\n');
                fprintf('  sSet     :: Data set i.e. EX, BY, PSI, etc\n');
                fprintf('  sSpecies :: Particle species i.e. EB, EP01, etc or for RAW, 01, 02, etc\n');
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
            sFile     = ['/',strrep(sFolder,'/','-'),'_',sTimeNExt,'.h5'];
            sLoad     = [obj.Path,'/',sFolder,sFile];
            
            % Check if datafile exists
            if iTime > iFiles
                fprintf(2, 'Error: Dump %d does not exist. Last dump is %d.\n', iTime, iFiles);
                return;
            end % if

            % Wraps reading in a try/catch because sometimes files are corrupt even if they exist.
            try
                if strcmpi(sType,'RAW')
                    stInfo  = h5info(sLoad,'/x1');
                    iSize   = stInfo.Dataspace.Size;
                    aReturn = zeros(iSize,6);

                    aReturn(:,1) = double(h5read(sLoad,'/x3'));
                    aReturn(:,2) = double(h5read(sLoad,'/x1'));
                    aReturn(:,3) = double(h5read(sLoad,'/x2'));
                    aReturn(:,4) = double(h5read(sLoad,'/p3'));
                    aReturn(:,5) = double(h5read(sLoad,'/p1'));
                    aReturn(:,6) = double(h5read(sLoad,'/p2'));
                else
                    aReturn = double(h5read(sLoad, ['/',sFolder]));
                end % if
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

                if strcmpi(sName,'RAW')
                    allDirs  = dir([obj.Path '/' sPath]);
                    allBeams = allDirs([allDirs(:).isdir]);
                    numBeams = numel(allBeams) - 2;
                    
                    stReturn.(sName).Info = struct('Path', '', 'Dirs', numBeams, 'Files', 0);

                    for b=1:numel(allBeams)
                        sBeam = allBeams(b).name;
                        if sum(ismember(sBeam, cExclude)) > 0
                            continue;
                        end % if

                        subFiles = dir([obj.Path '/' sPath '/' sBeam]);
                        numFiles = numel(subFiles) - 2;
                        sBName   = sprintf('EB%s',sBeam);
                        sBPath   = sprintf('%s/%s',sPath,sBeam);

                        stReturn.(sName).(sBName).Info = struct('Path', sBPath, 'Dirs', 0, 'Files', numFiles);
                    end % for
                else
                    subFiles = dir([obj.Path '/' sPath]);
                    numFiles = numel(subFiles) - 2;

                    stReturn.(sName).Info        = struct('Path', '', 'Dirs', 0, 'Files', 0);
                    stReturn.(sName).(sDir).Info = struct('Path', sPath, 'Dirs', 0, 'Files', numFiles);
                end % if
                
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
                    case 'R'
                        sDType    = 'RAW';
                        sDSet     = 'None';
                        sDSpecies = 'None';
                end % switch
                
                if ~isempty(sDType)
                    stSlice = fieldnames(obj.Elements.(sType));
                    for j=2:length(stSlice)
                        
                        if strcmpi(sDType,'RAW')
                            sDSpecies = stSlice{j};
                            sDSlice   = 'All';
                        else
                            sDSlice   = stSlice{j};
                        end % if

                        stData(iRow).Type    = sDType;
                        stData(iRow).Set     = sDSet;
                        stData(iRow).Species = sDSpecies;
                        stData(iRow).Slice   = sDSlice;
                        stData(iRow).Path    = obj.Elements.(sType).(stSlice{j}).Info.Path;
                        stData(iRow).Files   = obj.Elements.(sType).(stSlice{j}).Info.Files;
                        
                        stIndex.(sDType).(sDSet).(sDSpecies).(sDSlice) = iRow;
                        
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
