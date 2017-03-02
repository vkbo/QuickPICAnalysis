
%
%  Class Object :: Holds the QuickPIC config file
% ************************************************
%

classdef QPICConfig
    
    %
    % Properties
    %
    
    properties(GetAccess='public', SetAccess='public')

        Path       = '';     % Path to data directory
        File       = '';     % Config file within data directory
        Silent     = false;  % Set to true to disable command window output

    end % properties

    properties(GetAccess='public', SetAccess='private')
    
        Name       = '';     % Name of the loaded dataset
        Details    = {};     % Description of simulation. First lines of comments.
        Input      = {};     % Parsed input file
        Constants  = {};     % Constants
        Convert    = {};     % Unit conversion factors
        Simulation = {};     % Simulation variables
        EMFields   = {};     % Electro-magnetic field variables
        Particles  = {};     % Particle variables

    end % properties

    properties(GetAccess='public', SetAccess='private')

        Translate  = {};     % Container for Variables class
        NameLists  = {};     % Container for Fortran namelists
        
    end % properties

    %
    % Constructor
    %
    
    methods
        
        function obj = QPICConfig()
            
            % SI Constants
            obj.Constants.SI.SpeedOfLight       =  2.99792458e8;    % m/s (exact)
            obj.Constants.SI.ElectronMass       =  9.10938291e-31;  % kg
            obj.Constants.SI.ElementaryCharge   =  1.602176565e-19; % C
            obj.Constants.SI.VacuumPermitivity  =  8.854187817e-12; % F/m 
            obj.Constants.SI.VacuumPermeability =  1.2566370614e-6; % N/A^2
            obj.Constants.SI.Boltzmann          =  1.38064852e-23;  % J/K
            
            % Electron Volts
            obj.Constants.EV.ElectronMass       =  5.109989282e5;   % eV/c^2
            obj.Constants.EV.Boltzmann          =  8.6173324e-5;    % eV/K

            % CGS Constants
            obj.Constants.CGS.ElementaryCharge  =  4.80320425e-10;  % statC
            obj.Constants.CGS.Boltzmann         =  1.38064852e-16;  % erg/K
            
        end % function
        
    end % methods

    %
    % Setters an Getters
    %
    
    methods
        
        function obj = set.Path(obj, sPath)
            
            if ~isdir(sPath)
                return;
            end % if

            obj.Path = sPath;
            aPath    = strsplit(obj.Path, '/');
            obj.Name = aPath{end};
            sInput   = [sPath '/rpinput'];
            
            if exist(sInput,'file') == 2
                obj.File = 'rpinput';
            else
                fprintf(2, 'QPICConfig Error: Input file not found.\n');
            end % if

            if ~obj.Silent
                fprintf('Config file set: %s\n', obj.File);
            end % if

            obj = obj.fReadNameLists();
            obj = obj.fParseInputFile();

            %obj = obj.fGetSimulationVariables();
            %obj = obj.fGetEMFVariables();
            %obj = obj.fGetParticleVariables();

        end % function
                       
    end % methods
    
    %
    %  Config File Methods
    %
    
    methods(Access='private')
        
        function obj = fReadNameLists(obj)
            
            % Read file
            sFile  = fileread([obj.Path '/' obj.File]);
            
            % Get Simulation Description
            cLines = strsplit(sFile,'\n');
            iLines = numel(cLines);
            cDesc  = {};
            sClean = '';
            bStrip = false;

            for i=1:iLines
                sLine = strtrim(cLines{i});
                if ~isempty(sLine)
                    switch(sLine(1))
                        case '!'
                            cDesc{end+1} = strtrim(sLine(2:end));
                        case '-'
                            bStrip = ~bStrip;
                    end % switch
                end % if
                if isempty(sLine)
                    continue
                end % if
                if ~bStrip
                    if sLine(1) ~= '&' && sLine(1) ~= '/' && sLine(end) ~= ','
                        sLine = [sLine ','];
                    end % if
                    if sLine(1) ~= '-'
                        sClean = sprintf('%s%s\n',sClean,sLine);
                    end % if
                end % if
            end % for

            obj.Details = cDesc;
            
            % Clean-up
            sClean = regexprep(sClean,'\!.*?\n','\n');   % Removes ! style comments
            sClean = regexprep(sClean,'[\n|\r|\t]',' '); % Removes tabs and line breaks
            sClean = strrep(sClean,'.true.','1');        % Replace .true. with matlab logical true
            sClean = strrep(sClean,'.false.','0');       % Replace .false. with matlab logical false
            sClean = strrep(sClean,'"','''');            % Replace quote marks
            
            % Parse file and extract name lists

            sBuffer  = ' ';
            bQuote   = 0;
            bSection = 0;

            iNL = 1;
            stNameLists(iNL).Section  = [];
            stNameLists(iNL).NameList = [];

            for c=1:length(sClean)
                
                % Check if inside quote
                if sClean(c) == ''''
                    bQuote = ~bQuote;
                end % if

                % Add character to buffer, except spaces outside of quotes
                if ~(sClean(c) == ' ' && ~bQuote)
                    sBuffer = [sBuffer sClean(c)];
                end % if
                
                % If the character is a '&', the following word is a section until next space
                if sClean(c) == '&' && ~bQuote
                    bSection = 1;
                end % if
                if sClean(c) == ' ' && bSection
                    bSection = 0;
                    if length(sBuffer) > 1
                        sBuffer = strtrim(sBuffer);
                        stNameLists(iNL).Section = sBuffer(2:end);
                        sBuffer = ' '; % Reset buffer
                    end % if
                end % if
                
                % If the character is a '/', the buffer so far contains a name list
                if sClean(c) == '/' && ~bQuote
                    if length(sBuffer) > 1
                        sBuffer = strtrim(sBuffer(1:end-1));
                        if ~isempty(sBuffer)
                            if sBuffer(end) ~= ','
                                sBuffer = [sBuffer,','];
                            end % if
                        end % if
                        stNameLists(iNL).NameList = sBuffer;
                        sBuffer = ' '; % Reset buffer
                        iNL = iNL + 1; % Move to next section
                    end % if
                end % if
                
            end % for
            
            % Save the name list
            obj.NameLists = stNameLists;
            
        end % function

        function obj = fParseInputFile(obj)
            
            stInput = {};
            cGrp = {'simulation','beam','plasma','diagnostics'};
            iGrp = 1;
            
            aSim    = [];
            aBeam   = [];
            aPlasma = [];
            aDiag   = [];

            [~,iLMax] = size(obj.NameLists);
            for l=1:iLMax
                sSection = obj.NameLists(l).Section;
                switch(sSection)
                    case 'Pipeline'
                        iGrp = 1;
                    case 'Simulation_Sys'
                        iGrp = 1;
                    case 'Boundary'
                        iGrp = 1;
                    case 'Num_Beams'
                        iGrp = 1;
                    case 'Beam'
                        iGrp = 2;
                    case 'Plasma'
                        iGrp = 1;
                    case 'Species'
                        iGrp = 3;
                    case 'Simulation_time'
                        iGrp = 1;
                    case 'Potential_Diag'
                        iGrp = 4;
                    case 'Current_Diag'
                        iGrp = 4;
                    case 'Field_Diag'
                        iGrp = 4;
                    case 'Beam_Diag'
                        iGrp = 4;
                    case 'Plasma_Diag'
                        iGrp = 4;
                    case 'Beam_Phase_Space_Diag'
                        iGrp = 4;
                    case 'Restart_File'
                        iGrp = 1;
                    case 'debug'
                        iGrp = 1;
                end % switch
                switch(iGrp)
                    case 1; aSim    = [aSim l];
                    case 2; aBeam   = [aBeam l];
                    case 3; aPlasma = [aPlasma l];
                    case 4; aDiag   = [aDiag l];
                end % switch
            end % for

            % General Simulation Parameters
            for s=1:length(aSim)
                iSec = aSim(s);
                sSec = obj.NameLists(iSec).Section;
                sVal = obj.NameLists(iSec).NameList;
                stInput.(cGrp{1}).(sSec) = obj.fParseNameList(sVal);
            end % for
            
            % Beams
            for s=1:length(aBeam)
                iSec = aBeam(s);
                sSec = obj.NameLists(iSec).Section;
                sVal = sprintf('%s_%d',obj.NameLists(iSec).NameList,iSec);
                stInput.(cGrp{2}).(sSec) = obj.fParseNameList(sVal);
            end % for

            % Plasma
            for s=1:length(aPlasma)
                iSec = aPlasma(s);
                sSec = obj.NameLists(iSec).Section;
                sVal = sprintf('%s_%d',obj.NameLists(iSec).NameList,iSec);
                stInput.(cGrp{3}).(sSec) = obj.fParseNameList(sVal);
            end % for

            % Diagnostics
            for s=1:length(aDiag)
                iSec = aDiag(s);
                sSec = obj.NameLists(iSec).Section;
                sVal = obj.NameLists(iSec).NameList;
                stInput.(cGrp{4}).(sSec) = obj.fParseNameList(sVal);
            end % for
            % Save Struct
            obj.Input = stInput;
        
        end % function

    end % methods

    %
    %  Variables Methods
    %

    methods(Static, Access='private')
    
        function stReturn = fParseNameList(sData)
            
            stReturn = {};
            
            sBuffer = '';
            bQuote  = 0;
            bBrack  = 0;
            
            sVar = '';
            cVar = {};

            for c=1:length(sData)
                
                sBuffer = [sBuffer sData(c)];
                
                % Check if inside quote
                if sData(c) == ''''
                    bQuote = ~bQuote;
                end % if
                
                % Check if inside brackets
                if sData(c) == '('
                    bBrack = 1;
                end % if
                if sData(c) == ')'
                    bBrack = 0;
                end % if
                
                % After each variable follows an equal
                if sData(c) == '=' && ~bQuote && ~bBrack
                    if ~isempty(sVar)
                        try
                            evalc(['stReturn.',sVar,'={',strjoin(cVar,','),'}']);
                        catch
                            fprintf(2,'Error: Cannot parse variable ''%s'' in input file.\n',sVar);
                        end % try
                    end % if
                    sVar = sBuffer(1:end-1);
                    cVar = {};
                    sBuffer = '';
                end % if

                % After each value follows a comma
                if sData(c) == ',' && ~bQuote && ~bBrack
                    cVar{end+1} = sBuffer(1:end-1);
                    sBuffer = '';
                end % if

            end % for

            if ~isempty(sVar)
                try
                    evalc(['stReturn.',sVar,'={',strjoin(cVar,','),'}']);
                catch
                    fprintf(2,'Error: Cannot parse variable ''%s'' in input file.\n',sVar);
                end % try
            end % if
            
        end % function

    end % methods

end % classdef

