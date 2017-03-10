
%
%  Class Object :: Analyse Vecotrised Data Types
% ***********************************************
%  SubClass of QPICType
%
%  Description:
%    A class to analyse and handle QuickPIC data related to electric fields, magnetic fields, wakefields, and current.
%
%  Constructor:
%    oData   : QPICData object
%    sVector : What vector variable to load
%    Parameter pairs:
%      Units     : 'N', 'SI' or 'CGS'
%      X1Scale   : Unit scale on x1 axis. 'Auto', or specify metric unit
%      X2Scale   : Unit scale on x2 axis. 'Auto', or specify metric unit
%      X3Scale   : Unit scale on x3 axis. 'Auto', or specify metric unit
%      Scale     : Unit scale on all axes. 'Auto', or specify metric unit
%      Symmetric : Makes the X and Y axes symmetric around 0
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
%    Density  : Returns a dataset with a 2D matrix of the density of the field.
%    Lineout  : Returns a dataset with a 1D lineout of the density of the field.
%    Integral : Returns a dataset with the integrated field over an interval
%               of time dumps.
%

classdef QPICVector < QPICType

    %
    % Properties
    %

    properties(GetAccess='public', SetAccess='private')
        
        VectorVar   = '';  % Holds vector name
        VectorType  = '';  % Holds vector type
        VectorAxis  = 0;   % Holds vector axis
        VectorFac   = 1.0; % Vector scale factor
        VectorUnit  = 'N'; % Vector base unit
        VectorName  = '';  % Full name of the vector type
        VectorShort = '';  % Short name of the vector type
        VectorTex   = '';  % TeX name of the vector type
        
    end % properties

    %
    % Constructor
    %

    methods
        
        function obj = QPICVector(oData, sVector, varargin)
            
            % Call QPICType constructor
            obj@QPICType(oData, varargin{:});

            % Set Field
            sVector   = lower(sVector);
            sAxisName = '';
            cValid    = {'ex','ey','ez','bx','by','bz','wr','wz','jx','jy','jz'};
            if sum(ismember(sVector, cValid)) == 1
                obj.VectorVar  = sVector;
                obj.VectorType = sVector(1);
                switch sVector(2)
                    case 'x'
                        obj.VectorAxis = 2;
                        sAxisName      = 'Horizontal';
                    case 'y'
                        obj.VectorAxis = 3;
                        sAxisName      = 'Vertical';
                    case 'z'
                        obj.VectorAxis = 1;
                        sAxisName      = 'Longitudinal';
                    case 'r'
                        obj.VectorAxis = 4;
                        sAxisName      = 'Radial';
                end % switch
            else
                fprintf(2, 'Error: ''%s'' is not a recognised vector field. Using ''ez'' instead.\n', sVector);
                obj.VectorVar  = 'ez';
                obj.VectorType = 'e';
                obj.VectorAxis = 1;
            end % if
            
            if strcmpi(obj.Units,'SI')
                switch obj.VectorType
                    case 'e'
                        obj.VectorFac  = obj.Data.Config.Convert.SI.E0;
                        obj.VectorUnit = 'V/m';
                    case 'b'
                        obj.VectorFac  = obj.Data.Config.Convert.SI.B0;
                        obj.VectorUnit = 'T';
                    case 'w'
                        obj.VectorFac  = obj.Data.Config.Convert.SI.E0;
                        obj.VectorUnit = 'V/m';
                    case 'j'
                        obj.VectorFac  = obj.Data.Config.Convert.SI.JFac;
                        obj.VectorUnit = 'A';
                end % switch
            end % if
            
            % Names
            switch obj.VectorType
                case 'e'
                    obj.VectorName  = sprintf('%s E-field',sAxisName);
                    obj.VectorShort = sprintf('E%s',sVector(2));
                    obj.VectorTex   = sprintf('E_{%s}',sVector(2));
                case 'b'
                    obj.VectorName  = sprintf('%s B-field',sAxisName);
                    obj.VectorShort = sprintf('B%s',sVector(2));
                    obj.VectorTex   = sprintf('B_{%s}',sVector(2));
                case 'w'
                    obj.VectorName  = sprintf('%s Wakefield',sAxisName);
                    obj.VectorShort = sprintf('W%s',sVector(2));
                    obj.VectorTex   = sprintf('W_{%s}',sVector(2));
                case 'j'
                    obj.VectorName  = sprintf('%s Current',sAxisName);
                    obj.VectorShort = sprintf('J%s',sVector(2));
                    obj.VectorTex   = sprintf('J_{%s}',sVector(2));
            end % switch

        end % function
        
    end % methods

    %
    % Public Methods
    %
    
    methods(Access = 'public')
        
        function stReturn = Density2D(obj, sSlice)

            % Input/Output
            stReturn = {};
            
            if nargin < 2
                sSlice = '';
            end % if

            % Get Data and Parse it
            switch obj.VectorType
                case 'e'
                    aData = obj.Data.Data(obj.Time,'F',obj.VectorVar,'',sSlice);
                case 'b'
                    aData = obj.Data.Data(obj.Time,'F',obj.VectorVar,'',sSlice);
                case 'j'
                    aData = obj.Data.Data(obj.Time,'J',obj.VectorVar,'',sSlice);
                case 'w'
                    fprintf(2,'QPICVector.Density2D: For wakefields use Wakefield2D instead.\n');
                    return;
            end % switch
            if isempty(aData)
                fprintf(2,'QPICVector.Density2D: No data\n');
                return;
            end % if
            
            % If data is 3D, slice it
            if ndims(aData) == 3
                stData = obj.fParseGridData2D(aData);
                aData  = stData.Data;
                aHAxis = stData.HAxis;
                aVAxis = stData.VAxis;
                cAxes  = stData.Axes;
            else
                cAxes  = obj.fSliceOrientation(sSlice);
                aHAxis = obj.fGetBoxAxis(cAxes{1});
                aVAxis = obj.fGetBoxAxis(cAxes{2});
            end % if
            
            % Return Data
            stReturn.Data  = aData*obj.VectorFac;
            stReturn.HAxis = aHAxis;
            stReturn.VAxis = aVAxis;
            stReturn.Axes  = cAxes;
            stReturn.ZPos  = obj.fGetZPos();        
        
        end % function

        function stReturn = Lineout(obj, iStart, iAverage)

            % Input/Output
            stReturn = {};
            
            if nargin < 3
                iAverage = 1;
            end % if
            
            if nargin < 2
                iStart = 3;
            end % if
            
            % Get Data and Parse it
            aData = obj.Data.Data(obj.Time,'FLD',obj.FieldVar.Name,'');
            if isempty(aData)
                return;
            end % if

            stData = obj.fParseGridData1D(aData,iStart,iAverage);

            if isempty(stData)
                return;
            end % if

            % Return Data
            stReturn.Data   = stData.Data*obj.FieldFac;
            stReturn.HAxis  = stData.HAxis;
            stReturn.HRange = stData.HLim;
            stReturn.VRange = stData.VLim;
            stReturn.ZPos   = obj.fGetZPos();        
        
        end % function
        
        function stReturn = Wakefield2D(obj, sWakefield, varargin)
            
            % Input/Output
            stReturn = {};
            
            stWF = obj.Translate.Lookup(sWakefield,'Wakefield');
            if ~stWF.isWakefield
                fprintf(2, 'Error: ''%s'' is not a recognised wakefield. Using ''w1'' instead.\n', sWakefield);
                stWF = obj.Translate.Lookup('w1');
            end % if
            
            % Parse input
            oOpt = inputParser;
            addParameter(oOpt, 'GridDiag', {});
            parse(oOpt, varargin{:});
            stOpt = oOpt.Results;

            % Extract the beta values for the direction of simulation movement
            aMove = obj.Data.Config.Simulation.Moving;
            
            % Extract parallel e-field
            switch(stWF.Name)
                case 'w1'
                    aE    = obj.Data.Data(obj.Time,'FLD','e1','','GridDiag',stOpt.GridDiag);
                    bSign = false;
                case 'w2'
                    aE    = obj.Data.Data(obj.Time,'FLD','e2','','GridDiag',stOpt.GridDiag);
                    bSign = true;
                case 'w3'
                    aE    = obj.Data.Data(obj.Time,'FLD','e3','','GridDiag',stOpt.GridDiag);
                    bSign = true;
            end % switch
            
            % If parallel e-field doesn't exist, return empty
            if isempty(aE)
                return;
            end % if
            
            aB1 = 0.0;
            aB2 = 0.0;

            % Extract ortogonal b-fields – partial derivatives
            switch(stWF.Name)
                case 'w1'
                    if aMove(2) ~= 0.0
                        aB1 = aMove(2) * obj.Data.Data(obj.Time,'FLD','b3','','GridDiag',stOpt.GridDiag);
                    end % if
                    if aMove(3) ~= 0.0
                        aB2 = aMove(3) * obj.Data.Data(obj.Time,'FLD','b2','','GridDiag',stOpt.GridDiag);
                    end % if
                case 'w2'
                    if aMove(3) ~= 0.0
                        aB1 = aMove(3) * obj.Data.Data(obj.Time,'FLD','b1','','GridDiag',stOpt.GridDiag);
                    end % if
                    if aMove(1) ~= 0.0
                        aB2 = aMove(1) * obj.Data.Data(obj.Time,'FLD','b3','','GridDiag',stOpt.GridDiag);
                    end % if
                case 'w3'
                    if aMove(1) ~= 0.0
                        aB1 = aMove(1) * obj.Data.Data(obj.Time,'FLD','b2','','GridDiag',stOpt.GridDiag);
                    end % if
                    if aMove(2) ~= 0.0
                        aB2 = aMove(2) * obj.Data.Data(obj.Time,'FLD','b1','','GridDiag',stOpt.GridDiag);
                    end % if
            end % switch
            
            aW = aE + aB1 - aB2;
            
            % Slice the data
            stData = obj.fParseGridData2D(aW,'SignFlip',bSign,'GridDiag',stOpt.GridDiag);
            if isempty(stData)
                return;
            end % if
            
            dScale = obj.Data.Config.Convert.SI.E0;
            
            % Return Data
            stReturn.Data  = stData.Data*dScale;
            stReturn.HAxis = stData.HAxis;
            stReturn.VAxis = stData.VAxis;
            stReturn.Axes  = stData.Axes;
            stReturn.ZPos  = obj.fGetZPos();        
            
        end % function

        function stReturn = WFLineout(obj, sWakefield, iStart, iAverage)

            % Input/Output
            stReturn = {};
            
            stWF = obj.Translate.Lookup(sWakefield,'Wakefield');
            if ~stWF.isWakefield
                fprintf(2, 'Error: ''%s'' is not a recognised wakefield. Using ''w1'' instead.\n', sWakefield);
                stWF = obj.Translate.Lookup('w1');
            end % if
            
            % Extract the beta values for the direction of simulation movement
            aMove = obj.Data.Config.Simulation.Moving;
            
            % Extract parallel e-field
            switch(stWF.Name)
                case 'w1'
                    aE    = obj.Data.Data(obj.Time, 'FLD', 'e1', '');
                case 'w2'
                    aE    = obj.Data.Data(obj.Time, 'FLD', 'e2', '');
                case 'w3'
                    aE    = obj.Data.Data(obj.Time, 'FLD', 'e3', '');
            end % switch
            
            % If parallel e-field doesn't exist, return empty
            if isempty(aE)
                return;
            end % if
            
            aB1 = 0.0;
            aB2 = 0.0;

            % Extract ortogonal b-fields – partial derivatives
            switch(stWF.Name)
                case 'w1'
                    if aMove(2) ~= 0.0
                        aB1 = aMove(2) * obj.Data.Data(obj.Time, 'FLD', 'b3', '');
                    end % if
                    if aMove(3) ~= 0.0
                        aB2 = aMove(3) * obj.Data.Data(obj.Time, 'FLD', 'b2', '');
                    end % if
                case 'w2'
                    if aMove(3) ~= 0.0
                        aB1 = aMove(3) * obj.Data.Data(obj.Time, 'FLD', 'b1', '');
                    end % if
                    if aMove(1) ~= 0.0
                        aB2 = aMove(1) * obj.Data.Data(obj.Time, 'FLD', 'b3', '');
                    end % if
                case 'w3'
                    if aMove(1) ~= 0.0
                        aB1 = aMove(1) * obj.Data.Data(obj.Time, 'FLD', 'b2', '');
                    end % if
                    if aMove(2) ~= 0.o
                        aB2 = aMove(2) * obj.Data.Data(obj.Time, 'FLD', 'b1', '');
                    end % if
            end % switch
            
            aW = aE + aB1 - aB2;

            % Slice the data
            stData = obj.fParseGridData1D(aW,iStart,iAverage);
            if isempty(stData)
                return;
            end % if

            dScale = obj.Data.Config.Convert.SI.E0;

            % Return Data
            stReturn.Data   = stData.Data*dScale;
            stReturn.HAxis  = stData.HAxis;
            stReturn.HRange = stData.HLim;
            stReturn.VRange = stData.VLim;
            stReturn.ZPos   = obj.fGetZPos();        
        
        end % function
        
        function stReturn = Integral(obj, sStart, sStop, aRange)

            % Input/Output
            stReturn = {};

            if nargin < 2
                sStart = 'PStart';
            end % if

            if nargin < 3
                sStop = 'End';
            end % if

            if nargin < 4
                aRange = [];
            end % if

            iStart = obj.Data.StringToDump(sStart);
            iStop  = obj.Data.StringToDump(sStop);

            % Get simulation variables
            dTFac = obj.Data.Config.Convert.SI.TimeFac;
            dLFac = obj.Data.Config.Convert.SI.LengthFac;
            
            % Set axes
            aVAxis = [];
            aRAxis = [];
            aVLim  = [];
            sVUnit = 'N';
            sTUnit = 'm';

            switch(obj.FieldVar.Dim)

                case 1
                    dVFac  = obj.AxisFac(1);
                    sVUnit = obj.AxisUnits{1};
                    aVAxis = obj.fGetBoxAxis('x1');
                    aRAxis = obj.fGetBoxAxis('x2');
                    aVLim  = [fGetIndex(aVAxis, obj.X1Lim(1)*obj.AxisFac(1)) ...
                              fGetIndex(aVAxis, obj.X1Lim(2)*obj.AxisFac(1))];
                    
                    if isempty(aRange) || ~length(aRange) == 2
                        aRange = [3 3];
                    else
                        if aRange(1) < 1
                            aRange(1) = 1;
                        end % if
                        if aRange(2) > length(aRAxis)
                            aRange(2) = length(aRAxis);
                        end % if
                    end % if

                case 2
                    dVFac  = obj.AxisFac(2);
                    sVUnit = obj.AxisUnits{2};
                    aVAxis = obj.fGetBoxAxis('x2');
                    aRAxis = obj.fGetBoxAxis('x1');
                    aVLim  = [fGetIndex(aVAxis, obj.X2Lim(1)*obj.AxisFac(2)) ...
                              fGetIndex(aVAxis, obj.X2Lim(2)*obj.AxisFac(2))];

                    if isempty(aRange) || ~length(aRange) == 2
                        aRange = [1 10];
                    else
                        if aRange(1) < 1
                            aRange(1) = 1;
                        end % if
                        if aRange(2) > length(aRAxis)
                            aRange(2) = length(aRAxis);
                        end % if
                    end % if
                    
                case 3
                    return;

            end % switch
            
            aTAxis  = obj.fGetTimeAxis;
            aTAxis  = aTAxis(iStart+1:iStop+1);
            dTDiff  = aTAxis(end)-aTAxis(1);
            aVRange = [aVAxis(1) aVAxis(end)];
            aVAxis  = aVAxis(aVLim(1):aVLim(2));
            
            % Extract data
            aEnergy = zeros(length(aVAxis),length(aTAxis));
            for t=iStart:iStop
                
                aData = obj.Data.Data(t,'FLD',obj.FieldVar.Name,'');
                if isempty(aData)
                    return;
                end % if

                switch(obj.FieldVar.Name)
                    case 'e1'
                        aEnergy(:,t-iStart+1) = mean(aData(aVLim(1):aVLim(2),aRange(1):aRange(2)),2);
                    case 'e2'
                        aEnergy(:,t-iStart+1) = mean(aData(aRange(1):aRange(2),aVLim(1):aVLim(2)),1);
                end % switch
                
            end % for

            % Return data
            stReturn.Energy    = aEnergy*obj.FieldFac;
            stReturn.Integral  = cumtrapz(aEnergy,2)*obj.FieldFac*dTFac*dLFac;
            stReturn.GainFac   = 1/dTDiff;
            stReturn.VAxis     = aVAxis;
            stReturn.TAxis     = aTAxis;
            stReturn.VUnit     = sVUnit;
            stReturn.TUnit     = sTUnit;
            stReturn.AxisFac   = [1.0 dVFac];
            stReturn.AxisRange = [iStart iStop aVRange];
        
        end % function
        
    end % methods

end % classdef
