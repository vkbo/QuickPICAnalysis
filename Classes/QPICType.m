
%
%  Class Object :: SuperClass for QPICAnalysis Data Types
% ********************************************************
%

classdef QPICType

    %
    % Properties
    %
    
    properties(GetAccess='public', SetAccess='public')
        
        Time        = 1;                         % Current time (dumb number)
        X1Lim       = [];                        % Axes limits x1
        X2Lim       = [];                        % Axes limits x2
        X3Lim       = [];                        % Axes limits x3
        SliceAxis   = 3;                         % Slice axis (3D)
        Slice       = 0;                         % Slice coordinate (3D)
        Error       = false;                     % True if object failed to initialise

    end % properties

    properties(GetAccess='public', SetAccess='protected')
        
        Data        = {};                        % QPICData dataset
        Diagnostics = {};                        % Holds relevant diagnostics info
        Units       = 'N';                       % Units of axes
        AxisUnits   = {'N' 'N' 'N'};             % Units of axes
        AxisScale   = {'Auto' 'Auto' 'Auto'};    % Scale of axes
        AxisRange   = [0.0 0.0 0.0 0.0 0.0 0.0]; % Max and min of axes
        AxisFac     = [1.0 1.0 1.0];             % Axes scale factors
        Origin      = [0.0 0.0 0.0];             % Position defined as 0 on the grid axes
        SIOptions   = {};                        % Optional settings for SI Units
        Convert     = {};                        % Holds the struct of conversion factors

    end % properties
    
    properties(GetAccess='protected', SetAccess='protected')
        
        Translate   = {};                        % Lookup class for variables
        
    end % properties

    %
    % Constructor
    %

    methods
        
        function obj = QPICType(oData, varargin)
            
            % Set Data
            obj.Data = oData;

            % Read config
            dLFactor = obj.Data.Config.Convert.SI.LengthFac;
            aXMin    = obj.Data.Config.Simulation.XMin/dLFactor;
            aXMax    = obj.Data.Config.Simulation.XMax/dLFactor;

            % Read Input Parameters
            oOpt = inputParser;
            addParameter(oOpt, 'Units',     'N');
            addParameter(oOpt, 'Scale',     '');
            addParameter(oOpt, 'X1Scale',   'Auto');
            addParameter(oOpt, 'X2Scale',   'Auto');
            addParameter(oOpt, 'X3Scale',   'Auto');
            addParameter(oOpt, 'Symmetric', 'No');
            addParameter(oOpt, 'Tesla',     'No');
            addParameter(oOpt, 'JDensity',   'mm');
            addParameter(oOpt, 'QDensity',   'n0');
            parse(oOpt, varargin{:});
            stOpt = oOpt.Results;

            % Set Scale
            if isempty(stOpt.Scale)
                obj.AxisScale = {stOpt.X1Scale, stOpt.X2Scale, stOpt.X3Scale};
            else
                obj.AxisScale = {stOpt.Scale, stOpt.Scale, stOpt.Scale};
            end % if
            
            % Set Symmetric X/Y Axis
            if strcmpi(stOpt.Symmetric,'Yes')
                obj.Origin = [0.0 0.5*aXMax(2) 0.5*aXMax(3)];
                aXMin      = aXMin - obj.Origin;
                aXMax      = aXMax - obj.Origin;
            end % if
            
            % Units for B-fields
            if strcmpi(stOpt.Tesla,'Yes')
                obj.SIOptions.Tesla = true;
            else
                obj.SIOptions.Tesla = false;
            end % if
            
            % Metric Unit for Current Density (Squared)
            switch(lower(stOpt.JDensity))
                case 'm'
                    obj.SIOptions.JDensity = 1;
                case 'mm'
                    obj.SIOptions.JDensity = 2;
                case 'µm'
                    obj.SIOptions.JDensity = 3;
                case 'um'
                    obj.SIOptions.JDensity = 3;
                otherwise
                    obj.SIOptions.JDensity = 2;
            end % switch

            % Metric Unit for Charge Density (cubed)
            switch(lower(stOpt.QDensity))
                case 'n0'
                    obj.SIOptions.QDensity = 0;
                case 'm'
                    obj.SIOptions.QDensity = 1;
                case 'mm'
                    obj.SIOptions.QDensity = 2;
                case 'µm'
                    obj.SIOptions.QDensity = 3;
                case 'um'
                    obj.SIOptions.QDensity = 3;
                otherwise
                    obj.SIOptions.QDensity = 2;
            end % switch

            % Evaluate Units
            switch(lower(stOpt.Units))

                case 'si'
                    obj.Units = 'SI';
                    
                    [dX1Fac, sX1Unit]  = QPICTools.fLengthScale(obj.AxisScale{1}, 'm');
                    [dX2Fac, sX2Unit]  = QPICTools.fLengthScale(obj.AxisScale{2}, 'm');
                    [dX3Fac, sX3Unit]  = QPICTools.fLengthScale(obj.AxisScale{3}, 'm');
                    obj.AxisFac        = [dX1Fac dX2Fac dX3Fac]*dLFactor;
                    obj.AxisUnits      = {sX1Unit, sX2Unit, sX3Unit};
                    obj.AxisRange(1:2) = [aXMin(1) aXMax(1)]*obj.AxisFac(1);
                    obj.AxisRange(3:4) = [aXMin(2) aXMax(2)]*obj.AxisFac(2);
                    obj.AxisRange(5:6) = [aXMin(3) aXMax(3)]*obj.AxisFac(3);
                    obj.Convert        = obj.Data.Config.Convert.SI;

                case 'cgs'
                    obj.Units = 'CGS';
                    
                    [dX1Fac, sX1Unit]  = QPICTools.fLengthScale(obj.AxisScale{1}, 'cm');
                    [dX2Fac, sX2Unit]  = QPICTools.fLengthScale(obj.AxisScale{2}, 'cm');
                    [dX3Fac, sX3Unit]  = QPICTools.fLengthScale(obj.AxisScale{3}, 'cm');
                    obj.AxisFac        = [dX1Fac dX2Fac dX3Fac]*dLFactor;
                    obj.AxisUnits      = {sX1Unit, sX2Unit, sX3Unit};
                    obj.AxisRange(1:2) = [aXMin(1) aXMax(1)]*obj.AxisFac(1);
                    obj.AxisRange(3:4) = [aXMin(2) aXMax(2)]*obj.AxisFac(2);
                    obj.AxisRange(5:6) = [aXMin(3) aXMax(3)]*obj.AxisFac(3);
                    obj.Convert        = obj.Data.Config.Convert.CGS;

                otherwise
                    obj.Units = 'N';

                    obj.AxisFac        = [1.0, 1.0, 1.0];
                    obj.AxisUnits      = {'c/ω', 'c/ω', 'c/ω'};
                    obj.AxisRange      = [aXMin(1) aXMax(1) aXMin(2) aXMax(2) aXMin(3) aXMax(3)];
                    obj.Convert        = obj.Data.Config.Convert.Norm;

            end % switch

            % Set defult axis limits
            obj.X1Lim = [aXMin(1) aXMax(1)]*obj.AxisFac(1);
            obj.X2Lim = [aXMin(2) aXMax(2)]*obj.AxisFac(2);
            obj.X3Lim = [aXMin(3) aXMax(3)]*obj.AxisFac(3);
            
            % Set default slice for 3D
            obj.SliceAxis = 3;
            obj.Slice     = 0.0;

        end % function

    end % methods

    %
    % Setters and Getters
    %

    methods
        
        function obj = set.Time(obj, iTime)
            
            dTMax = obj.Data.Config.Simulation.TMax;
            dDT   = obj.Data.Config.Simulation.TimeStep;
            iTMax = floor(dTMax/dDT);

            if iTime < 1
                fprintf(2,'QPICType: Invalid time %d. Minimum time is 1.\n',iTime);
                iTime = 1;
            end % if
            
            if iTime > iTMax
                fprintf(2,'QPICType: Invalid time %d. Maxmum time is %d.\n',iTime,iTMax);
                iTime = iTMax;
            end % if
            
            obj.Time = iTime;
            
        end % function
        
        function obj = set.X1Lim(obj, aLim)

            dXMin = obj.AxisRange(1);
            dXMax = obj.AxisRange(2);

            if length(aLim) ~= 2
                fprintf(2, 'Error: x1 limit needs to be a vector of dimension 2.\n');
                return;
            end % if

            if aLim(2) < aLim(1)
                fprintf(2, 'Error: second value must be larger than first value.\n');
                return;
            end % if

            if aLim(1) < dXMin || aLim(1) > dXMax || aLim(2) < dXMin || aLim(2) > dXMax
                fprintf('Warning: Lim input is out of range. Range is %.2f–%.2f %s.\n',dXMin,dXMax,obj.AxisUnits{1});
                aLim = obj.AxisRange(1:2);
            end % if

            obj.X1Lim = aLim/obj.AxisFac(1);

        end % function
         
        function obj = set.X2Lim(obj, aLim)
 
            dXMin = obj.AxisRange(3);
            dXMax = obj.AxisRange(4);

            if length(aLim) ~= 2
                fprintf(2, 'Error: x2 limit needs to be a vector of dimension 2.\n');
                return;
            end % if

            if aLim(2) < aLim(1)
                fprintf(2, 'Error: second value must be larger than first value.\n');
                return;
            end % if
            
            if aLim(1) < dXMin || aLim(1) > dXMax || aLim(2) < dXMin || aLim(2) > dXMax
                fprintf('Warning: Lim input is out of range. Range is %.2f–%.2f %s.\n',dXMin,dXMax,obj.AxisUnits{2});
                aLim = obj.AxisRange(3:4);
            end % if

            obj.X2Lim = aLim/obj.AxisFac(2);
             
        end % function
 
        function obj = set.X3Lim(obj, aLim)

            dXMin = obj.AxisRange(5);
            dXMax = obj.AxisRange(6);

            if length(aLim) ~= 2
                fprintf(2, 'Error: x3 limit needs to be a vector of dimension 2.\n');
                return;
            end % if

            if aLim(2) < aLim(1)
                fprintf(2, 'Error: second value must be larger than first value.\n');
                return;
            end % if

            if aLim(1) < dXMin || aLim(1) > dXMax || aLim(2) < dXMin || aLim(2) > dXMax
                fprintf('Warning: Lim input is out of range. Range is %.2f–%.2f %s.\n',dXMin,dXMax,obj.AxisUnits{3});
                aLim = obj.AxisRange(5:6);
            end % if

            obj.X3Lim = aLim/obj.AxisFac(3);

        end % function
        
        function obj = set.SliceAxis(obj, iAxis)
            
            if isempty(iAxis)
                return;
            end % if
            
            iAxis = floor(iAxis);
            
            if iAxis > 0 && iAxis <= 3
                obj.SliceAxis = iAxis;
            else
                fprintf(2, 'Error: Not a proper axis.\n');
            end % if
            
            obj.Slice = 0.0;
            
        end % function

        function obj = set.Slice(obj, dSlice)
            
            if isempty(dSlice)
                return;
            end % if

            aAxis     = obj.fGetBoxAxis(sprintf('x%d',obj.SliceAxis));
            obj.Slice = QPICTools.fGetIndex(aAxis,dSlice);
            
        end % function

    end % methods

    %
    % Public Methods
    %
    
    methods(Access='public')

        function sReturn = PlasmaPosition(obj)

            dTFac   = obj.Data.Config.Convert.SI.TimeFac;
            dLFac   = obj.Data.Config.Convert.SI.LengthFac;
            dSimPos = obj.Time*dTFac*dLFac;

            [dTemp,sUnit] = QPICTools.fAutoScale(dSimPos,'m',1e-6);
            dScale  = dTemp/dSimPos;

            sReturn = sprintf('z=%0.2f %s',dSimPos*dScale,sUnit);

        end % function

        function sReturn = SlicePosition(obj,sSlice)

            dLFac  = obj.Data.Config.Convert.SI.LengthFac;
            aXMax  = obj.Data.Config.Simulation.XMax;
            aGrid  = obj.Data.Config.Simulation.Grid;
            aDelta = aXMax ./ aGrid;

            [sSlice,iDim,iOrth] = QPICTools.fCheckSlice(sSlice);

            switch(iDim)
                case 2
                    sSlice = QPICTools.fLabelAxis(iOrth,true);
                    dSlice = obj.Diagnostics.Slices(iOrth) + 0.5*aDelta(iOrth) - obj.Origin(iOrth)*dLFac;
                case 3
                    sSlice = QPICTools.fLabelAxis(obj.SliceAxis,true);
                    dSlice = (obj.Slice-0.5)*aDelta(obj.SliceAxis) - obj.Origin(obj.SliceAxis)*dLFac;
            end % switch
            
            [dTemp,sUnit] = QPICTools.fAutoScale(dSlice,'m',1e-6);
            dScale  = dTemp/dSlice;

            sReturn = sprintf('%s=%0.2f %s',sSlice,dSlice*dScale,sUnit);

        end % function

    end % methods

    %
    % Protected Methods
    %
    
    methods(Access='protected')
        
        function stReturn = fParseGridData1D(obj, aData, iStart, iAverage)

            % Input/Output
            stReturn = {};
            
            iDim       = ndims(aData);
            iSliceAxis = obj.SliceAxis;
            
            switch iSliceAxis
                case 1
                    sHAxis = 'x2';
                    sVAxis = 'x3';
                    aHAxis = obj.fGetBoxAxis('x2');
                    aVAxis = obj.fGetBoxAxis('x3');
                    aHLim  = [obj.X2Lim(1)*obj.AxisFac(2), obj.X2Lim(2)*obj.AxisFac(2)];
                    aVLim  = [obj.X3Lim(1)*obj.AxisFac(3), obj.X3Lim(2)*obj.AxisFac(3)];
                case 2
                    sHAxis = 'x1';
                    sVAxis = 'x3';
                    aHAxis = obj.fGetBoxAxis('x1');
                    aVAxis = obj.fGetBoxAxis('x3');
                    aHLim  = [obj.X1Lim(1)*obj.AxisFac(1), obj.X1Lim(2)*obj.AxisFac(1)];
                    aVLim  = [obj.X3Lim(1)*obj.AxisFac(3), obj.X3Lim(2)*obj.AxisFac(3)];
                case 3
                    sHAxis = 'x1';
                    sVAxis = 'x2';
                    aHAxis = obj.fGetBoxAxis('x1');
                    aVAxis = obj.fGetBoxAxis('x2');
                    aHLim  = [obj.X1Lim(1)*obj.AxisFac(1), obj.X1Lim(2)*obj.AxisFac(1)];
                    aVLim  = [obj.X2Lim(1)*obj.AxisFac(2), obj.X2Lim(2)*obj.AxisFac(2)];
            end % switch

            if iDim == 3
                switch iSliceAxis
                    case 1
                        aData  = squeeze(aData(obj.Slice,:,:));
                    case 2
                        aData  = squeeze(aData(:,obj.Slice,:));
                    case 3
                        aData  = squeeze(aData(:,:,obj.Slice));
                end % switch
            end % if
            
            % Get H-Limits
            iHMin  = fGetIndex(aHAxis, aHLim(1));
            iHMax  = fGetIndex(aHAxis, aHLim(2));
            
            % Get V-Limits
            iVN    = numel(aVAxis);

            % Crop Dataset
            iEnd   = iStart+iAverage-1;
            aData  = squeeze(mean(aData(iHMin:iHMax,iStart:iEnd),2));
            aHAxis = aHAxis(iHMin:iHMax);
            
            % Return Data
            stReturn.Data  = aData;
            stReturn.HAxis = aHAxis;
            stReturn.HLim  = [iHMin iHMax];
            stReturn.VLim  = [aVAxis(iStart) aVAxis(iEnd+1)];
            stReturn.Axes  = {sHAxis,sVAxis};

        end % function
        
        function stReturn = fParseGridData2D(obj, aData)
            
            % Input/Output
            stReturn = {};
            
            iDim       = ndims(aData);
            iSliceAxis = obj.SliceAxis;

            if iDim == 2
                return;
            end % if

            switch iSliceAxis
                case 1
                    sHAxis = 'x2';
                    sVAxis = 'x3';
                    aHAxis = obj.fGetBoxAxis('x2');
                    aVAxis = obj.fGetBoxAxis('x3');
                    aHLim  = [obj.X2Lim(1)*obj.AxisFac(2), obj.X2Lim(2)*obj.AxisFac(2)];
                    aVLim  = [obj.X3Lim(1)*obj.AxisFac(3), obj.X3Lim(2)*obj.AxisFac(3)];
                case 2
                    sHAxis = 'x1';
                    sVAxis = 'x3';
                    aHAxis = obj.fGetBoxAxis('x1');
                    aVAxis = obj.fGetBoxAxis('x3');
                    aHLim  = [obj.X1Lim(1)*obj.AxisFac(1), obj.X1Lim(2)*obj.AxisFac(1)];
                    aVLim  = [obj.X3Lim(1)*obj.AxisFac(3), obj.X3Lim(2)*obj.AxisFac(3)];
                case 3
                    sHAxis = 'x1';
                    sVAxis = 'x2';
                    aHAxis = obj.fGetBoxAxis('x1');
                    aVAxis = obj.fGetBoxAxis('x2');
                    aHLim  = [obj.X1Lim(1)*obj.AxisFac(1), obj.X1Lim(2)*obj.AxisFac(1)];
                    aVLim  = [obj.X2Lim(1)*obj.AxisFac(2), obj.X2Lim(2)*obj.AxisFac(2)];
            end % switch

            if iDim == 3
                switch iSliceAxis
                    case 1
                        aData  = squeeze(aData(:,:,obj.Slice));
                    case 2
                        aData  = squeeze(aData(obj.Slice,:,:));
                    case 3
                        aData  = squeeze(aData(:,obj.Slice,:));
                end % switch
            end % if

            % Get Limits
            iHMin = QPICTools.fGetIndex(aHAxis, aHLim(1));
            iHMax = QPICTools.fGetIndex(aHAxis, aHLim(2));
            iVMin = QPICTools.fGetIndex(aVAxis, aVLim(1));
            iVMax = QPICTools.fGetIndex(aVAxis, aVLim(2));

            % Crop Dataset
            aData  = aData(iVMin:iVMax,iHMin:iHMax);
            aHAxis = aHAxis(iHMin:iHMax);
            aVAxis = aVAxis(iVMin:iVMax);
            
            % Return Data
            stReturn.Data  = aData;
            stReturn.HAxis = aHAxis;
            stReturn.VAxis = aVAxis;
            stReturn.HLim  = [iHMin iHMax];
            stReturn.VLim  = [iVMin iVMax];
            stReturn.Axes  = {sHAxis,sVAxis};

        end % function
        
        function aReturn = fGetTimeAxis(obj)
            
            iDumps  = obj.Data.SimData.MaxFiles;
            
            dTFac   = obj.Convert.TimeFac;
            dLFac   = obj.Convert.LengthFac;
            
            aReturn = linspace(0.0, dTFac*iDumps, iDumps)*dLFac;
            
        end % function

        function aReturn = fGetBoxAxis(obj, sAxis)
            
            dLFac     = obj.Convert.LengthFac;
            [~,iAxis] = QPICTools.fTranslateAxis(sAxis);
            
            if iAxis > 0 && iAxis < 4
                dXMin  = obj.Data.Config.Simulation.XMin(iAxis)/dLFac - obj.Origin(iAxis);
                dXMax  = obj.Data.Config.Simulation.XMax(iAxis)/dLFac - obj.Origin(iAxis);
                iNX    = obj.Data.Config.Simulation.Grid(iAxis);
                dScale = obj.AxisFac(iAxis);
            else
                dXMin  = 0.0;
                dXMax  = 0.0;
                iNX    = 0;
                dScale = 1.0;
            end % if

            aReturn = linspace(dXMin, dXMax, iNX)*dScale;
            
        end % function
        
        function dReturn = fGetZPos(obj)
            
            dLFactor = obj.Convert.LengthFac;
            dTFactor = obj.Convert.TimeFac;
            dReturn  = obj.Time*dTFactor*dLFactor;
            
        end % function
        
        function dReturn = fGetXPos(obj,sAxis,iIndex)
            
            aAxis   = obj.fGetBoxAxis(sAxis);
            dReturn = aAxis(iIndex);
            
        end % function
        
        function bReturn = fError(obj)
            
            bReturn = false;
            
            if obj.Error
                
                fprintf(2, 'QPICType object not initialised properly.\n');
                bReturn = true;
                
            end % if
            
        end % function

        function aReturn = fMomentumToEnergy(obj, aMomentum, dMass)
            
            dEMass  = obj.Data.Config.Constants.EV.ElectronMass;
            aReturn = sqrt(aMomentum.^2 + 1)*dMass*dEMass;
            
        end % function
        
        function aRaw = fRawToXi(obj, aRaw)

            dOffset = obj.Data.Config.Simulation.TimeStep*obj.Data.Config.Simulation.NDump*obj.Time;
            aRaw(:,1) = aRaw(:,1) - dOffset;
            
        end % function

    end % methods

    %
    %  Static Methods
    %

    methods(Static, Access='protected')
     
        function aReturn = fPruneRaw(aRaw, iVar, dMin, dMax)
            
            % By default do nothing.
            aReturn = aRaw;
            
            if iVar < 1 || iVar > 6 || isempty(dMin) && isempty(dMax)
                return;
            end % if
            
            if ~isempty(dMin)
                aInd = aReturn(:,iVar) < dMin;
                aReturn(aInd,:) = [];
            end % if

            if ~isempty(dMax)
                aInd = aReturn(:,iVar) > dMax;
                aReturn(aInd,:) = [];
            end % if
            
        end % function
        
        function bReturn = fValidSlice(sSlice)
            
            switch lower(sSlice)
                case ''
                    bReturn = true;
                case 'xy'
                    bReturn = true;
                case 'xz'
                    bReturn = true;
                case 'yz'
                    bReturn = true;
                otherwise
                    bReturn = false;
            end % switch
            
        end % function

        function cReturn = fSliceOrientation(sSlice)
            
            switch lower(sSlice)
                case 'xy'
                    cReturn = {'x2' 'x3'};
                case 'xz'
                    cReturn = {'x1' 'x2'};
                case 'yz'
                    cReturn = {'x1' 'x3'};
                otherwise
                    cReturn = {'' ''};
            end % switch
            
        end % function
        
    end % methods
    
end % classdef
