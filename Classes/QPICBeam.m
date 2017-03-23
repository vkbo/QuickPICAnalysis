
%
%  Class Object :: QPICBeam
% **************************
%  SubClass of QPICType
%
%  Description:
%    A class to analyse and handle QuickPIC data related to the beam particles.
%
%  Constructor:
%    oData : QPICData object
%    iBeam : What beam to load
%    Parameter pairs:
%      Units   : 'N', 'SI' or 'CGS'
%      X1Scale : Unit scale on x1 axis. 'Auto', or specify metric unit
%      X2Scale : Unit scale on x2 axis. 'Auto', or specify metric unit
%      X3Scale : Unit scale on x3 axis. 'Auto', or specify metric unit
%      Scale   : Unit scale on all axes. 'Auto', or specify metric unit
%
%  Set Methods:
%    Time      : Set time dump for dataset. Default is 0.
%    X1Lim     : 2D array of limits for x1 axis. Default is full box.
%    X2Lim     : 2D array of limits for x2 axis. Default is full box.
%    X3Lim     : 2D array of limits for x3 axis. Default is full box.
%    SliceAxis : 2D slice axis for 3D data
%    Slice     : 2D slice coordinate for 3D data
%
%  Public Methods:
%    PhaseSpace : Returns a dataset with information on beam emittance
%                 or emittance related properties.
%

classdef QPICBeam < QPICType

    %
    % Properties
    %

    properties(GetAccess='public', SetAccess='public')
        
        BeamVar   = 1;   % Internal beam name
        BeamNum   = 1;   % Beam number
        BeamFac   = 1.0; % Sample factor of the beam dummp
        BeamName  = '';  % Full name of the beam
        BeamShort = '';  % Short name of the beam
        BeamTex   = '';  % TeX name of the beam
        BeamConf  = {};  % Beam simulation parameters

    end % properties

    %
    % Constructor
    %

    methods
        
        function obj = QPICBeam(oData, iBeam, varargin)
            
            % Call OsirisType constructor
            obj@QPICType(oData, varargin{:});
            
            if iBeam < 1 || iBeam > 99
                fprintf(2,'QPICPhaseSpace: Invalid beam number. Must be between 1 and 99.\n');
                iBeam = 1;
            end % if
            
            obj.Diagnostics = obj.Data.Config.Diag.RAW;
            obj.BeamNum     = iBeam;
            obj.BeamVar     = sprintf('EB%02d',iBeam);
            obj.BeamConf    = obj.Data.Config.Beam.(sprintf('Beam%02d',iBeam));
            obj.BeamFac     = double(obj.Diagnostics.Sample);
            
            
            % Guess the beam type
            dCharge       = obj.BeamConf.Charge;
            dMass         = obj.BeamConf.Mass;

            obj.BeamName  = sprintf('Beam %1d',iBeam);
            obj.BeamShort = sprintf('EB%1d',iBeam);
            obj.BeamTex   = sprintf('e^{%1d}',iBeam);
            switch dCharge
                case -1
                    if dMass == 1.0
                        obj.BeamName  = 'Electron Beam';
                        obj.BeamShort = 'EB';
                        obj.BeamTex   = 'e^{-}';
                    end % if
                case 1
                    if dMass == 1.0
                        obj.BeamName  = 'Positron Beam';
                        obj.BeamShort = 'EB+';
                        obj.BeamTex   = 'e^{+}';
                    end % if
                    if dMass >= 1830.0 && dMass <= 1840.0
                        obj.BeamName  = 'Proton Beam';
                        obj.BeamShort = 'PB';
                        obj.BeamTex   = 'p^{+}';
                    end % if
            end % switch

        end % function
        
    end % methods

    %
    % Public Methods
    %
    
    methods(Access='public')
        
        function stReturn = PhaseSpace(obj, varargin)
            %
            %  Function: QPICPhaseSpace.PhaseSpace
            % *************************************
            %
            %  Options:
            % ==========
            %  Histogram [Yes/No]
            %          Whether to use fixed bins like in a histogram. Default behaviour is to deposit the charge on
            %          a sample grid. Default 'no'.
            %  Grid [2-vector, integer]
            %          The dimension of the hisorgram or deposit grid.
            %  Dimension [Rad/RadToX/X/Y]
            %          Which dimension to calculate the angle on.
            %          [Rad]    The radial axis alone.
            %          [RadToX] Projecting r and theta onto x by sampling theta.
            %          [X], [Y] The x or y axis alone.
            %          Defaults to RadToX for cylindrical and X for cartesian simulations.
            %  Slice [2-vector, float]
            %          Excludes all macro particles outside the defind limits. Units defined by class setting Units, and
            %          Scale parameter. Default is none, i.e. all particles.
            %  SliceAxis [1:6]
            %          Which axis to slice along. Setting is ignored if Slice is not defined. Defaults to 1.
            %
            %  Outputs:
            % ==========
            %  Raw        :: The macroparticles used in the calculation
            %  X          :: The x axis values
            %  XUnit      :: The unit of the x axis
            %  XPrime     :: The x prime axis (angle)
            %  XPrimeUnit :: The unit of the x prime axis
            %  Charge     :: The charge of the macroparticles used
            %  Covariance :: The covariance matrix
            %  ERMS       :: The RMS emittance
            %  ENorm      :: The normalised emittance
            %  BeamGamma  :: The mean gamma of the particles
            %  Alpha      :: The alpha twiss parameter
            %  Beta       :: The beta twiss parameter
            %  Gamma      :: The gamma twiss parameter
            %

            % Input/Output
            stReturn       = {};
            stReturn.Error = '';

            % Check that the object is initialised
            if obj.fError
                return;
            end % if

            oOpt = inputParser;
            addParameter(oOpt, 'Histogram', 'No');
            addParameter(oOpt, 'Grid',      [1000 1000]);
            addParameter(oOpt, 'Dimension', '');
            addParameter(oOpt, 'Slice',     '');
            addParameter(oOpt, 'SliceAxis', 1);
            parse(oOpt, varargin{:});
            stOpt = oOpt.Results;
            
            switch lower(stOpt.Dimension)
                case 'x'
                    iDim = 2;
                case 'y'
                    iDim = 3;
                otherwise
                    iDim = 2;
            end % switch

            aRaw = obj.Data.Data(obj.Time,'RAW','',obj.BeamVar,'');
            if isempty(aRaw)
                stReturn.Error = 'No initial data.';
                return;
            end % if
            
            % Slice data if requested
            if ~isempty(stOpt.Slice)
                if stOpt.SliceAxis >= 1 && stOpt.SliceAxis <= 6 && numel(stOpt.Slice) == 2
                    if stOpt.SliceAxis <= 3
                        aSlice = stOpt.Slice/obj.AxisFac(stOpt.SliceAxis);
                    end % if
                    if stOpt.SliceAxis > 3 && stOpt.SliceAxis <= 6
                        dFac   = obj.Data.Config.Constants.EV.ElectronMass;
                        dFac   = dFac*obj.BeamConf.Mass;
                        aSlice = stOpt.Slice/dFac;
                    end % if
                    aRaw = obj.fPruneRaw(aRaw,stOpt.SliceAxis,aSlice(1),aSlice(2));
                end % if
            end % if

            if isempty(aRaw)
                stReturn.Error = 'No data after cut.';
                return;
            end % if

            % Calculate Emittance

            aPz     = aRaw(:,4);
            aPx     = aRaw(:,iDim+3);
            aX      = aRaw(:,iDim)*obj.AxisFac(iDim);
            aX      = aX-mean(aX);

            aXP     = tan(aPx./aPz)*1e3;
            aCov    = cov(aX, aXP);
            dGamma  = mean(aPz);
            dERMS   = sqrt(det(aCov));
            dENorm  = sqrt(det(aCov))*dGamma;

            % Return
            
            stReturn.Raw        = aRaw;
            stReturn.X          = aX;
            stReturn.XUnit      = obj.AxisUnits{1};
            stReturn.XPrime     = aXP;
            stReturn.XPrimeUnit = 'mrad';
            stReturn.Covariance = aCov;
            stReturn.ERMS       = dERMS;
            stReturn.ENorm      = dENorm;
            stReturn.BeamGamma  = dGamma;
            stReturn.ZPos       = obj.fGetZPos();
            stReturn.Count      = length(aX);
            stReturn.RawCount   = length(aRaw);
            stReturn.Charge     = length(aX)*obj.BeamFac*obj.BeamConf.SimCharge;

            % Twiss parameters
            stReturn.Alpha      = aCov(1,2)/dERMS;
            stReturn.Beta       = aCov(1,1)/dERMS;
            stReturn.Gamma      = aCov(2,2)/dERMS;

            if strcmpi(stOpt.Histogram, 'No')
                return;
            end % if
            
            aHist   = zeros(stOpt.Grid(1),stOpt.Grid(2));
            
            dXMax   = max(abs(aX));
            dXPMax  = max(abs(aXP));
            dXMin   = -dXMax;
            dXPMin  = -dXPMax;
            dQ      = obj.BeamFac*obj.BeamConf.SimCharge;
            
            % Check if angle is 0 (no emittance)
            if dXPMax == 0

                stReturn.Error = 'Maximum particle angle is 0.';
                stReturn.Hist  = [];
                stReturn.HAxis = [];
                stReturn.VAxis = [];
                stReturn.Count = 0;
                
                return;

            end % if
                
            dDX  = dXMax/((stOpt.Grid(1)-2)/2);
            dDXP = dXPMax/((stOpt.Grid(2)-2)/2);
            aM   = floor(aX/dDX)+(stOpt.Grid(1)/2)+1;
            aN   = floor(aXP/dDXP)+(stOpt.Grid(2)/2)+1;

            for i=1:length(aM)
                aHist(aM(i),aN(i)) = aHist(aM(i),aN(i)) + dQ;
            end % for

            stReturn.Hist  = abs(transpose(aHist))*1e9;
            stReturn.HAxis = linspace(dXMin,dXMax,stOpt.Grid(1));
            stReturn.VAxis = linspace(dXPMin,dXPMax,stOpt.Grid(2));
            
        end % function

        function stReturn = SlicedPhaseSpace(obj, varargin)

            % Input/Output
            stReturn       = {};
            stReturn.Error = '';

            % Check that the object is initialised
            if obj.fError
                return;
            end % if

            oOpt = inputParser;
            addParameter(oOpt, 'Dimension', 'x');
            addParameter(oOpt, 'Lim',       []);
            addParameter(oOpt, 'Smooth',    4.0);
            addParameter(oOpt, 'MinStat',   100);
            parse(oOpt, varargin{:});
            stOpt = oOpt.Results;
            
            switch lower(stOpt.Dimension)
                case 'x'
                    iDim = 2;
                case 'y'
                    iDim = 3;
                otherwise
                    iDim = 2;
            end % switch

            aRaw = obj.Data.Data(obj.Time,'RAW','',obj.BeamVar,'');
            if isempty(aRaw)
                stReturn.Error = 'No initial data.';
                return;
            end % if
        
            dLFac = obj.Convert.LengthFac;
            dZMax = obj.Data.Config.Simulation.XMax(1);
            iGrid = obj.Data.Config.Simulation.Grid(1);
            dDz   = dZMax/double(iGrid)/dLFac;
            
            if isempty(stOpt.Lim)
                iMinZ = 1;
                iMaxZ = iGrid;
            else
                iMinZ = stOpt.Lim(1);
                iMaxZ = stOpt.Lim(2);
            end % if
            
            iNz = iMaxZ-iMinZ+1;
            
            aEmN = zeros(iNz,1);
            aEmG = zeros(iNz,1);
            aEx  = zeros(iNz,1);
            aNP  = zeros(iNz,1);

            for s=1:iNz

                iZ = s+iMinZ-1;

                dLMin = (iZ - stOpt.Smooth*0.5)*dDz;
                dLMax = (iZ + stOpt.Smooth*0.5)*dDz;

                aRawT = aRaw;
                aRawT((aRawT(:,1) < dLMin | aRawT(:,1) > dLMax),:) = [];
                aNP(s) = numel(aRawT(:,1));

                if aNP(s) < stOpt.MinStat
                    aEx(s) = 1;
                    continue;
                end % if

                aPz  = aRawT(:,4);
                aPx  = aRawT(:,iDim+3);
                aX   = aRawT(:,iDim)*obj.AxisFac(iDim);
                aX   = aX-mean(aX);

                aXP  = tan(aPx./aPz)*1e3;
                aCov = cov(aX, aXP);
                dMPz = mean(aPz);
                dEm  = real(sqrt(det(aCov)));

                aEmG(s) = dEm;
                aEmN(s) = dEm*dMPz;

            end % for
            
            aAxis = obj.fGetBoxAxis('x1');
            aAxis = aAxis(iMinZ:iMaxZ);
            
            stReturn.ERMS       = aEmG;
            stReturn.ENorm      = aEmN;
            stReturn.Axis       = aAxis;
            stReturn.XUnit      = obj.AxisUnits{1};
            stReturn.XPrimeUnit = 'mrad';
            stReturn.Count      = aNP;
            stReturn.Excluded   = aEx;

        end % function
    
    end % methods

end % classdef
