
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
        Beam       = {};     % Beam settings, but no diagnostics
        Plasma     = {};     % Plasma settings, but no diagnostics
        Diag       = {};     % All simulation diagnostics

    end % properties

    properties(GetAccess='private', SetAccess='private')

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
    % Setters and Getters
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

            obj = obj.fGetSimulationVariables();
            obj = obj.fGetParticleVariables();
            obj = obj.fGetDiagVariables();

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
                sSection = lower(obj.NameLists(l).Section);
                switch(sSection)
                    case 'pipeline'
                        iGrp = 1;
                    case 'simulation_sys'
                        iGrp = 1;
                    case 'boundary'
                        iGrp = 1;
                    case 'num_beams'
                        iGrp = 1;
                    case 'beam'
                        iGrp = 2;
                    case 'plasma'
                        iGrp = 1;
                    case 'species'
                        iGrp = 3;
                    case 'simulation_time'
                        iGrp = 1;
                    case 'potential_diag'
                        iGrp = 4;
                    case 'current_diag'
                        iGrp = 4;
                    case 'field_diag'
                        iGrp = 4;
                    case 'beam_diag'
                        iGrp = 4;
                    case 'plasma_diag'
                        iGrp = 4;
                    case 'beam_phase_space_diag'
                        iGrp = 4;
                    case 'restart_file'
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
                sSec = lower(obj.NameLists(iSec).Section);
                sVal = lower(obj.NameLists(iSec).NameList);
                stInput.(cGrp{1}).(sSec) = obj.fParseNameList(sVal);
            end % for
            
            % Beams
            for s=1:length(aBeam)
                iSec = aBeam(s);
                sSec = sprintf('%s_%d',lower(obj.NameLists(iSec).Section),s);
                sVal = obj.NameLists(iSec).NameList;
                stInput.(cGrp{2}).(sSec) = obj.fParseNameList(sVal);
            end % for

            % Plasma
            for s=1:length(aPlasma)
                iSec = aPlasma(s);
                sSec = sprintf('%s_%d',lower(obj.NameLists(iSec).Section),s);
                sVal = obj.NameLists(iSec).NameList;
                stInput.(cGrp{3}).(sSec) = obj.fParseNameList(sVal);
            end % for

            % Diagnostics
            for s=1:length(aDiag)
                iSec = aDiag(s);
                sSec = lower(obj.NameLists(iSec).Section);
                sVal = lower(obj.NameLists(iSec).NameList);
                stInput.(cGrp{4}).(sSec) = obj.fParseNameList(sVal);
            end % for
            % Save Struct
            obj.Input = stInput;
        
        end % function

    end % methods

    %
    %  Static Methods
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
                    sVar = lower(sBuffer(1:end-1));
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

        function aReturn = fArrayPad(aData, aTemplate)
            
            aTemplate(1:length(aData)) = aData;
            aReturn = aTemplate;
            
        end % function

    end % methods

    %
    %  Variable Methods
    %

    methods(Access = 'private')

        function obj = fGetSimulationVariables(obj)
            
            % Constants
            dC          = obj.Constants.SI.SpeedOfLight;
            dECharge    = obj.Constants.SI.ElementaryCharge;
            dEChargeCGS = obj.Constants.CGS.ElementaryCharge;
            dEMass      = obj.Constants.SI.ElectronMass;
            dEpsilon0   = obj.Constants.SI.VacuumPermitivity;
            dMu0        = obj.Constants.SI.VacuumPermeability;

            %
            % Main Simulation Variables
            %
            
            % Plasma Density
            try
                dN0 = double(obj.Input.simulation.plasma.plasma_density{1})*1.0e6;
            catch
                dN0 = 1.0;
            end % try
            
            % Plasma Frequency and Wavelength
            dOmegaP  = sqrt((dN0 * dECharge^2) / (dEMass * dEpsilon0));
            dLambdaP = 2*pi * dC / dOmegaP;
            
            % Coordinates
            try
                indX  = obj.Input.simulation.simulation_sys.indx{1};
                indY  = obj.Input.simulation.simulation_sys.indy{1};
                indZ  = obj.Input.simulation.simulation_sys.indz{1};
                aGrid = [2^indZ 2^indX 2^indY];
            catch
                aGrid = [1 1 1];
            end % try
            iDim  = 3;
            aGrid = obj.fArrayPad(aGrid, [0 0 0]);
            
            % Grid
            sCoords      = 'cartesian';
            bCylindrical = false;

            % Time Step
            try
                dTimeStep = double(obj.Input.simulation.simulation_time.dt{1});
            catch
                dTimeStep = 0.0;
            end % try

            % Start Time
            dTMin = 0.0;

            % End Time
            try
                dTMax = double(obj.Input.simulation.simulation_time.tend{1});
            catch
                dTMax = 0.0;
            end % try

            % Box Size
            aXMin = [0.0 0.0 0.0];
            try
                dMaxX = double(obj.Input.simulation.simulation_sys.box_x{1});
                dMaxY = double(obj.Input.simulation.simulation_sys.box_y{1});
                dMaxZ = double(obj.Input.simulation.simulation_sys.box_z{1});
                aXMax = [dMaxZ dMaxX dMaxY]*1e-6;
            catch
                aXMax = 0.0;
            end % try
            aXMax = obj.fArrayPad(aXMax, [0.0 0.0 0.0]);
            
            % Moving
            aMove = [1 0 0];

            % Save Plasma Variables for Simulation
            obj.Simulation.N0          = dN0;
            obj.Simulation.OmegaP      = dOmegaP;
            obj.Simulation.LambdaP     = dLambdaP;
            
            % Save Other Variables
            obj.Simulation.Coordinates = sCoords;
            obj.Simulation.Cylindrical = bCylindrical;
            obj.Simulation.Dimensions  = iDim;
            obj.Simulation.Grid        = aGrid;
            obj.Simulation.NDump       = 1;
            
            % Save Time Variables
            obj.Simulation.TimeStep    = dTimeStep;
            obj.Simulation.TMin        = dTMin;
            obj.Simulation.TMax        = dTMax;

            % Save Space Variables
            obj.Simulation.XMin        = aXMin;
            obj.Simulation.XMax        = aXMax;
            obj.Simulation.Moving      = aMove;

            
            %
            % Conversion Factors
            %
            
            % Geometry
            dLFactor = dC / dOmegaP;

            % Electric and Magnetic Field
            dSIE0 = dEMass * dC^3 * dOmegaP * dMu0*dEpsilon0 / dECharge;
            dSIB0 = dEMass * dC^2 * dOmegaP * dMu0*dEpsilon0 / dECharge;

            % Save Conversion Variables
            obj.Convert.SI.E0        = dSIE0;
            obj.Convert.SI.B0        = dSIB0;
            obj.Convert.SI.LengthFac = dLFactor;
            obj.Convert.SI.TimeFac   = dTimeStep;
            

            %
            % Particle and Charge Conversion
            %
            
            dDX1  = (aXMax(1) - aXMin(1))/aGrid(1);
            dDX2  = (aXMax(2) - aXMin(2))/aGrid(2);
            dDX3  = (aXMax(3) - aXMin(3))/aGrid(3);
            
            dPFac = dN0;          % Density is relative to N0
            dPFac = dPFac*dDX1;   % Scale for z cell size
            dPFac = dPFac*dDX2;   % Scale for x cell size
            dPFac = dPFac*dDX3;   % Scale for y cell size

            obj.Convert.Norm.ChargeFac   = 1.0;                % Already in units of n0
            obj.Convert.Norm.ParticleFac = dPFac/dLFactor^3;   % Particles per cell in units of c/omega_p
            obj.Convert.SI.ChargeFac     = dPFac*dECharge;     % Charge per m^2
            obj.Convert.SI.ParticleFac   = dPFac;              % Particles per m^2
            obj.Convert.CGS.ChargeFac    = dPFac*dEChargeCGS;  % Charge per cm^2
            obj.Convert.CGS.ParticleFac  = dPFac*1e-4;         % Particles per cm^2

            
            %
            % Current Conversion
            %

            aJFac    = [1.0 1.0 1.0]; % In normalised units

            aJFacSI  = aJFac     * dPFac*dECharge*dC;
            aJFacSI  = aJFacSI  ./ [dDX1 dDX2 dDX3];
            aJFacCGS = aJFac     * dPFac*dEChargeCGS*dC;
            aJFacCGS = aJFacCGS ./ [dDX1 dDX2 dDX3];

            obj.Convert.Norm.JFac = aJFac;
            obj.Convert.SI.JFac   = aJFacSI;
            obj.Convert.CGS.JFac  = aJFacCGS;

        end % function

        function obj = fGetParticleVariables(obj)
            
            % Get Number of Beams
            try
                iBeams = int64(obj.Input.simulation.num_beams.nbeams{1});
            catch
                iBeams = 0;
            end % try

            % Get Number of Plasmas
            try
                iPlasmas = int64(obj.Input.simulation.plasma.nspecies{1});
            catch
                iPlasmas = 0;
            end % try
            
            % Save Particle Values
            obj.Beam.NBeams     = iBeams;
            obj.Plasma.NPlasmas = iPlasmas;
            
            %
            % Loop All Beams
            %
            
            for p=1:iBeams
                
                sRawName = sprintf('beam_%d',p);
                sNewName = sprintf('Beam%02d',p);

                % Beam Evolution
                try
                    iBeamEvol = int64(obj.Input.beam.(sRawName).beam_evolution{1});
                catch
                    iBeamEvol = 1;
                end % try

                % Beam Particles
                try
                    iNPX   = int64(obj.Input.beam.(sRawName).npx{1});
                    iNPY   = int64(obj.Input.beam.(sRawName).npy{1});
                    iNPZ   = int64(obj.Input.beam.(sRawName).npz{1});
                    aNPart = [iNPZ iNPX iNPY];
                catch
                    aNPart = [1 1 1];
                end % try
                
                % Beam Charge
                try
                    dCharge = double(obj.Input.beam.(sRawName).charge{1});
                catch
                    dCharge = 1.0;
                end % try

                % Beam Mass
                try
                    dMass = double(obj.Input.beam.(sRawName).mass{1});
                catch
                    dMass = 1.0;
                end % try

                % Beam Gamma
                try
                    dGamma = double(obj.Input.beam.(sRawName).gamma{1});
                catch
                    dGamma = 1.0;
                end % try

                % Beam Particles
                try
                    iCount = int64(obj.Input.beam.(sRawName).num_particle{1});
                catch
                    iCount = 1;
                end % try

                % Beam Parameters
                try
                    aParam = cell2mat(obj.Input.beam.(sRawName).parameter_array);
                    aParam = padarray(aParam, [5 3], 0.0, 'post');
                    aParam = aParam(1:5,1:3);
                catch
                    aParam = zeros(5,3);
                end % try
                
                % Calculate Additional Variables
                dECharge    = obj.Constants.SI.ElementaryCharge;
                iSimTotal   = prod(aNPart);
                dBeamCharge = dCharge * dECharge * double(iCount);
                dSimCharge  = dBeamCharge/double(iSimTotal);

                % Save Variables
                obj.Beam.(sNewName).BeamEvol     = iBeamEvol;
                obj.Beam.(sNewName).SimParticles = aNPart;
                obj.Beam.(sNewName).SimTotal     = iSimTotal;
                obj.Beam.(sNewName).SimCharge    = dSimCharge;
                obj.Beam.(sNewName).Charge       = dCharge;
                obj.Beam.(sNewName).Mass         = dMass;
                obj.Beam.(sNewName).Gamma        = dGamma;
                obj.Beam.(sNewName).NumParticles = iCount;
                obj.Beam.(sNewName).BeamCharge   = dBeamCharge;
                obj.Beam.(sNewName).Position     = aParam(1,[3 1 2])*1e-6;
                obj.Beam.(sNewName).Profile      = aParam(2,[3 1 2])*1e-6;
                obj.Beam.(sNewName).Emittance    = aParam(3,[3 1 2])*1e-6;
                obj.Beam.(sNewName).CentroidX    = aParam(4,:)*1e-6;
                obj.Beam.(sNewName).CentroidY    = aParam(5,:)*1e-6;

            end % for

            %
            % Loop All Plasmas
            %
            
            for p=1:iPlasmas
                
                sRawName = sprintf('species_%d',p);
                sNewName = sprintf('Plasma%02d',p);

                % Beam Particles
                try
                    iNP2 = int64(obj.Input.plasma.(sRawName).np2{1});
                catch
                    iNP2 = 1;
                end % try
                
                % Beam Charge
                try
                    dCharge = double(obj.Input.plasma.(sRawName).charge{1});
                catch
                    dCharge = 1.0;
                end % try

                % Beam Mass
                try
                    dMass = double(obj.Input.plasma.(sRawName).mass{1});
                catch
                    dMass = 1.0;
                end % try

                % Save Variables
                obj.Plasma.(sNewName).SliceParticles = iNP2;
                obj.Plasma.(sNewName).Charge         = dCharge;
                obj.Plasma.(sNewName).Mass           = dMass;

            end % for

        end % function
        
        function obj = fGetDiagVariables(obj)
            
            cDF = {};
            cDS = {};
            
            %
            % E-Field Diagnostics
            %
            
            % Dump
            try
                iDump      = int64(obj.Input.diagnostics.field_diag.dfe{1});
                cDF{end+1} = 'EField';
            catch
                iDump      = 0;
            end % try
            
            % Dump Slices
            try
                iDumpSlice = int64(obj.Input.diagnostics.field_diag.dfeslice{1});
                dX0        = double(obj.Input.diagnostics.field_diag.ex0{1});
                dY0        = double(obj.Input.diagnostics.field_diag.ey0{1});
                dZ0        = double(obj.Input.diagnostics.field_diag.ez0{1});
                aSlices    = [dZ0 dX0 dY0];
                cDS{end+1} = 'EField';
            catch
                iDumpSlice = 0;
                aSlices    = [0 0 0];
            end % try
            
            % Save Variables
            obj.Diag.EField.Dump      = iDump;
            obj.Diag.EField.DumpSlice = iDumpSlice;
            obj.Diag.EField.Slices    = aSlices*1e-6;

            %
            % B-Field Diagnostics
            %

            % Dump
            try
                iDump      = int64(obj.Input.diagnostics.field_diag.dfb{1});
                cDF{end+1} = 'BField';
            catch
                iDump      = 0;
            end % try
            
            % Dump Slices
            try
                iDumpSlice = int64(obj.Input.diagnostics.field_diag.dfeslice{1});
                dX0        = double(obj.Input.diagnostics.field_diag.bx0{1});
                dY0        = double(obj.Input.diagnostics.field_diag.by0{1});
                dZ0        = double(obj.Input.diagnostics.field_diag.bz0{1});
                aSlices    = [dZ0 dX0 dY0];
                cDS{end+1} = 'BField';
            catch
                iDumpSlice = 0;
                aSlices    = [0 0 0];
            end % try
            
            % Save Variables
            obj.Diag.BField.Dump      = iDump;
            obj.Diag.BField.DumpSlice = iDumpSlice;
            obj.Diag.BField.Slices    = aSlices*1e-6;

            %
            % Beam Diagnostics
            %
            
            % Dump
            try
                iDump      = int64(obj.Input.diagnostics.beam_diag.dfqeb{1});
                cDF{end+1} = 'Beam';
            catch
                iDump      = 0;
            end % try
            
            % Slices
            try
                iDumpSlice = int64(obj.Input.diagnostics.beam_diag.dfqebslice{1});
                dX0        = double(obj.Input.diagnostics.beam_diag.qebx0{1});
                dY0        = double(obj.Input.diagnostics.beam_diag.qeby0{1});
                dZ0        = double(obj.Input.diagnostics.beam_diag.qebz0{1});
                aSlices    = [dZ0 dX0 dY0];
                cDS{end+1} = 'Beam';
            catch
                iDumpSlice = 0;
                aSlices    = [0 0 0];
            end % try
            
            % Save Variables
            obj.Diag.Beam.Dump      = iDump;
            obj.Diag.Beam.DumpSlice = iDumpSlice;
            obj.Diag.Beam.Slices    = aSlices*1e-6;

            %
            % Plasma Diagnostics
            %
            
            % Dump
            try
                iDump      = int64(obj.Input.diagnostics.plasma_diag.dfqep{1});
                cDF{end+1} = 'Plasma';
            catch
                iDump      = 0;
            end % try
            
            % Slices
            try
                iDumpSlice = int64(obj.Input.diagnostics.plasma_diag.dfqepslice{1});
                dX0        = double(obj.Input.diagnostics.plasma_diag.qepx0{1});
                dY0        = double(obj.Input.diagnostics.plasma_diag.qepy0{1});
                dZ0        = double(obj.Input.diagnostics.plasma_diag.qepz0{1});
                aSlices    = [dZ0 dX0 dY0];
                cDS{end+1} = 'Plasma';
            catch
                iDumpSlice = 0;
                aSlices    = [0 0 0];
            end % try
            
            % Save Variables
            obj.Diag.Plasma.Dump      = iDump;
            obj.Diag.Plasma.DumpSlice = iDumpSlice;
            obj.Diag.Plasma.Slices    = aSlices*1e-6;

            %
            % Current Diagnostics
            %
            
            % Dump
            try
                iDump      = int64(obj.Input.diagnostics.current_diag.dfjp{1});
                cDF{end+1} = 'Current';
            catch
                iDump      = 0;
            end % try
            
            % Slices
            try
                iDumpSlice = int64(obj.Input.diagnostics.current_diag.dfjpslice{1});
                dX0        = double(obj.Input.diagnostics.current_diag.jpx0{1});
                dY0        = double(obj.Input.diagnostics.current_diag.jpy0{1});
                dZ0        = double(obj.Input.diagnostics.current_diag.jpz0{1});
                aSlices    = [dZ0 dX0 dY0];
                cDS{end+1} = 'Current';
            catch
                iDumpSlice = 0;
                aSlices    = [0 0 0];
            end % try
            
            % Save Variables
            obj.Diag.Current.Dump      = iDump;
            obj.Diag.Current.DumpSlice = iDumpSlice;
            obj.Diag.Current.Slices    = aSlices*1e-6;

            %
            % Potential Diagnostics
            %
            
            % Dump
            try
                iDump      = int64(obj.Input.diagnostics.potential_diag.dfpsi{1});
                cDF{end+1} = 'Potential';
            catch
                iDump      = 0;
            end % try
            
            % Slices
            try
                iDumpSlice = int64(obj.Input.diagnostics.potential_diag.dfpsislice{1});
                dX0        = double(obj.Input.diagnostics.potential_diag.psix0{1});
                dY0        = double(obj.Input.diagnostics.potential_diag.psiy0{1});
                dZ0        = double(obj.Input.diagnostics.potential_diag.psiz0{1});
                aSlices    = [dZ0 dX0 dY0];
                cDS{end+1} = 'Potential';
            catch
                iDumpSlice = 0;
                aSlices    = [0 0 0];
            end % try
            
            % Save Variables
            obj.Diag.Potential.Dump      = iDump;
            obj.Diag.Potential.DumpSlice = iDumpSlice;
            obj.Diag.Potential.Slices    = aSlices*1e-6;
            
            %
            %  Diagnostics Lists
            %
            
            obj.Diag.Diagnostics = cDF;
            obj.Diag.Slices      = cDS;
            
            %
            % RAW Dump Diagnostics
            %
            
            % Dump Enabled
            try
                bDump = int64(obj.Input.diagnostics.beam_phase_space_diag.dump_pha_beam{1});
            catch
                bDump = 0;
            end % try
            
            % Dump Time Step
            try
                iDump = int64(obj.Input.diagnostics.beam_phase_space_diag.dfpha_beam{1});
            catch
                iDump = 0;
            end % try
            
            % Dump Sample
            try
                iSample = int64(obj.Input.diagnostics.beam_phase_space_diag.dsample_beam{1});
            catch
                iSample = 0;
            end % try
            
            % Save Variables
            obj.Diag.RAW.Enabled  = bDump;
            obj.Diag.RAW.TimeStep = iDump;
            obj.Diag.RAW.Sample   = iSample;
            
        end % function

    end % methods

end % classdef

