
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
        ParticleFac = 1.0;                       % Q-to-particles factor
        ChargeFac   = 1.0;                       % Q-to-charge factor
        BoxOffset   = 0.0;                       % Start of the box in simulation
        XOrigin     = 0.0;                       % Position defined as 0 on the x axis
        YOrigin     = 0.0;                       % Position defined as 0 on the y axis
        SIOptions   = {};                        % Optional settings for SI Units

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
                obj.XOrigin = 0.5*aXMax(2);
                obj.YOrigin = 0.5*aXMax(3);
                aXMin(2:3)  = aXMin(2:3) - [obj.XOrigin obj.YOrigin];
                aXMax(2:3)  = aXMax(2:3) - [obj.XOrigin obj.YOrigin];
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
                    
                    [dX1Fac, sX1Unit]  = obj.fLengthScale(obj.AxisScale{1}, 'm');
                    [dX2Fac, sX2Unit]  = obj.fLengthScale(obj.AxisScale{2}, 'm');
                    [dX3Fac, sX3Unit]  = obj.fLengthScale(obj.AxisScale{3}, 'm');
                    obj.AxisFac        = [dX1Fac dX2Fac dX3Fac]*dLFactor;
                    obj.AxisUnits      = {sX1Unit, sX2Unit, sX3Unit};
                    obj.AxisRange(1:2) = [aXMin(1) aXMax(1)]*obj.AxisFac(1);
                    obj.AxisRange(3:4) = [aXMin(2) aXMax(2)]*obj.AxisFac(2);
                    obj.AxisRange(5:6) = [aXMin(3) aXMax(3)]*obj.AxisFac(3);
                    
                    obj.ParticleFac    = obj.Data.Config.Convert.SI.ParticleFac;
                    obj.ChargeFac      = obj.Data.Config.Convert.SI.ChargeFac;

                case 'cgs'
                    obj.Units = 'CGS';
                    
                    [dX1Fac, sX1Unit]  = obj.fLengthScale(obj.AxisScale{1}, 'cm');
                    [dX2Fac, sX2Unit]  = obj.fLengthScale(obj.AxisScale{2}, 'cm');
                    [dX3Fac, sX3Unit]  = obj.fLengthScale(obj.AxisScale{3}, 'cm');
                    obj.AxisFac        = [dX1Fac dX2Fac dX3Fac]*dLFactor;
                    obj.AxisUnits      = {sX1Unit, sX2Unit, sX3Unit};
                    obj.AxisRange(1:2) = [aXMin(1) aXMax(1)]*obj.AxisFac(1);
                    obj.AxisRange(3:4) = [aXMin(2) aXMax(2)]*obj.AxisFac(2);
                    obj.AxisRange(5:6) = [aXMin(3) aXMax(3)]*obj.AxisFac(3);
                    
                    obj.ParticleFac    = obj.Data.Config.Convert.CGS.ParticleFac;
                    obj.ChargeFac      = obj.Data.Config.Convert.CGS.ChargeFac;

                otherwise
                    obj.Units = 'N';

                    obj.AxisFac     = [1.0, 1.0, 1.0];
                    obj.AxisUnits   = {'c/ω', 'c/ω', 'c/ω'};
                    obj.AxisRange   = [aXMin(1) aXMax(1) aXMin(2) aXMax(2) aXMin(3) aXMax(3)];

                    obj.ParticleFac = obj.Data.Config.Convert.Norm.ParticleFac;
                    obj.ChargeFac   = obj.Data.Config.Convert.Norm.ChargeFac;

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
            
            if iTime < 1
                iTime = 1;
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
            sSlice = QPICTools.fCheckSlice(sSlice);
            aXMax  = obj.Data.Config.Simulation.XMax;
            aGrid  = obj.Data.Config.Simulation.Grid;
            aDelta = aXMax ./ aGrid;
            
            switch(sSlice)
                case 'XY'
                    sSlice = '\xi';
                    dSlice = obj.Diagnostics.Slices(1) + 0.5*aDelta(1);
                case 'YZ'
                    sSlice = 'x';
                    dSlice = obj.Diagnostics.Slices(2) + 0.5*aDelta(2) - obj.XOrigin*dLFac;
                case 'XZ'
                    sSlice = 'y';
                    dSlice = obj.Diagnostics.Slices(3) + 0.5*aDelta(3) - obj.XOrigin*dLFac;
                otherwise
                    % 3D data
                    switch(obj.SliceAxis)
                        case 1
                            sSlice = '\xi';
                            dSlice = (obj.Slice-0.5)*aDelta(1);
                        case 2
                            sSlice = 'x';
                            dSlice = (obj.Slice-0.5)*aDelta(2) - obj.XOrigin*dLFac;
                        case 3
                            sSlice = 'y';
                            dSlice = (obj.Slice-0.5)*aDelta(3) - obj.YOrigin*dLFac;
                    end % switch
            end % switch
            
            [dTemp,sUnit] = QPICTools.fAutoScale(dSlice,'m',1e-6);
            dScale  = dTemp/dSlice;

            sReturn = sprintf('%s=%0.2f %s',sSlice,dSlice*dScale,sUnit);

        end % function

    end % methods

    %
    % Private Methods
    %
    
    methods(Access='private')
        
        function [dScale, sUnit] = fLengthScale(~, sToUnit, sFromUnit)

            dScale = 1.0;
            sUnit  = 'm';

            if nargin < 2
                sFromUnit = 'm';
            end % if

            switch(lower(sFromUnit))
                case 'pm'
                    dScale = dScale * 1.0e-12;
                case 'å'
                    dScale = dScale * 1.0e-10;
                case 'nm'
                    dScale = dScale * 1.0e-9;
                case 'um'
                    dScale = dScale * 1.0e-6;
                case 'µm'
                    dScale = dScale * 1.0e-6;
                case 'mm'
                    dScale = dScale * 1.0e-3;
                case 'cm'
                    dScale = dScale * 1.0e-2;
                case 'm'
                    dScale = dScale * 1.0;
                case 'km'
                    dScale = dScale * 1.0e3;
            end % switch

            switch(lower(sToUnit))
                case 'pm'
                    dScale = dScale * 1.0e12;
                    sUnit  = 'pm';
                case 'å'
                    dScale = dScale * 1.0e10;
                    sUnit  = 'Å';
                case 'nm'
                    dScale = dScale * 1.0e9;
                    sUnit  = 'nm';
                case 'um'
                    dScale = dScale * 1.0e6;
                    sUnit  = 'µm';
                case 'µm'
                    dScale = dScale * 1.0e6;
                    sUnit  = 'µm';
                case 'mm'
                    dScale = dScale * 1.0e3;
                    sUnit  = 'mm';
                case 'cm'
                    dScale = dScale * 1.0e2;
                    sUnit  = 'cm';
                case 'm'
                    dScale = dScale * 1.0;
                    sUnit  = 'm';
                case 'km'
                    dScale = dScale * 1.0e-3;
                    sUnit  = 'km';
            end % switch

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
            
            dTFac   = obj.Data.Config.Convert.SI.TimeFac;
            dLFac   = obj.Data.Config.Convert.SI.LengthFac;
            
            aReturn = linspace(0.0, dTFac*iDumps, iDumps)*dLFac;
            
        end % function

        function aReturn = fGetBoxAxis(obj, sAxis)
            
            dLFac = obj.Data.Config.Convert.SI.LengthFac;
            sAxis = obj.fTranslateAxis(sAxis);

            switch sAxis
                case 'x1'
                    dXMin  = obj.Data.Config.Simulation.XMin(1)/dLFac;
                    dXMax  = obj.Data.Config.Simulation.XMax(1)/dLFac;
                    iNX    = obj.Data.Config.Simulation.Grid(1);
                    dScale = obj.AxisFac(1);
                case 'x2'
                    dXMin  = obj.Data.Config.Simulation.XMin(2)/dLFac - obj.XOrigin;
                    dXMax  = obj.Data.Config.Simulation.XMax(2)/dLFac - obj.XOrigin;
                    iNX    = obj.Data.Config.Simulation.Grid(2);
                    dScale = obj.AxisFac(2);
                case 'x3'
                    dXMin  = obj.Data.Config.Simulation.XMin(3)/dLFac - obj.YOrigin;
                    dXMax  = obj.Data.Config.Simulation.XMax(3)/dLFac - obj.YOrigin;
                    iNX    = obj.Data.Config.Simulation.Grid(3);
                    dScale = obj.AxisFac(3);
                otherwise
                    dXMin  = 0.0;
                    dXMax  = 0.0;
                    iNX    = 0;
                    dScale = 1.0;
            end % switch

            aReturn = linspace(dXMin, dXMax, iNX)*dScale;
            
        end % function
        
        function dReturn = fGetZPos(obj)
            
            dLFactor = obj.Data.Config.Convert.SI.LengthFac;
            dTFactor = obj.Data.Config.Convert.SI.TimeFac;
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
        
        function sReturn = fTranslateAxis(sAxis)
            
            switch lower(sAxis)
                case 'x'
                    sReturn = 'x2';
                case 'y'
                    sReturn = 'x3';
                case 'z'
                    sReturn = 'x1';
                otherwise
                    sReturn = sAxis;
            end % switch
            
        end % function
        
        function sReturn = fLabelAxis(sAxis)
            
            switch lower(sAxis)
                case 'x'
                    sReturn = 'x';
                case 'y'
                    sReturn = 'y';
                case 'z'
                    sReturn = '\xi';
                case 'r'
                    sReturn = 'r';
                case 'x1'
                    sReturn = '\xi';
                case 'x2'
                    sReturn = 'x';
                case 'x3'
                    sReturn = 'y';
                case 't'
                    sReturn = 'z';
                otherwise
                    sReturn = sAxis;
            end % switch
            
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
