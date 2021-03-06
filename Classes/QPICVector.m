
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
            
            % Values
            dC = obj.Data.Config.Constants.SI.SpeedOfLight;

            % Set Field
            sVector   = lower(sVector);
            sAxisName = '';
            cValid    = {'ex','ey','ez','bx','by','bz','wx','wy','wz','wr','jpx','jpy','jpz'};
            if sum(ismember(sVector, cValid)) == 1
                obj.VectorVar  = sVector;
                obj.VectorType = sVector(1);
                switch sVector(end)
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
                        obj.VectorFac   = obj.Data.Config.Convert.SI.E0;
                        obj.VectorUnit  = 'V/m';
                        obj.Diagnostics = obj.Data.Config.Diag.EField;
                    case 'b'
                        if obj.SIOptions.Tesla
                            obj.VectorFac  = obj.Data.Config.Convert.SI.B0;
                            obj.VectorUnit = 'T';
                        else
                            obj.VectorFac  = obj.Data.Config.Convert.SI.B0*dC;
                            obj.VectorUnit = 'V/mc';
                        end % if
                        obj.Diagnostics = obj.Data.Config.Diag.BField;
                    case 'w'
                        obj.VectorFac   = obj.Data.Config.Convert.SI.E0;
                        obj.VectorUnit  = 'V/m';
                        obj.Diagnostics = obj.Data.Config.Diag.EField;
                    case 'j'
                        switch(obj.SIOptions.JDensity)
                            case 1
                                obj.VectorFac  = obj.Data.Config.Convert.SI.JFac(obj.VectorAxis);
                                obj.VectorUnit = 'A/m^2';
                            case 2
                                obj.VectorFac  = obj.Data.Config.Convert.SI.JFac(obj.VectorAxis)*1e6;
                                obj.VectorUnit = 'A/mm^2';
                            case 3
                                obj.VectorFac  = obj.Data.Config.Convert.SI.JFac(obj.VectorAxis)*1e12;
                                obj.VectorUnit = 'A/µm^2';
                        end % switch
                        obj.Diagnostics = obj.Data.Config.Diag.Current;
                end % switch
            end % if
            
            % Names
            switch obj.VectorType
                case 'e'
                    obj.VectorName  = sprintf('%s E-Field',sAxisName);
                    obj.VectorShort = sprintf('E%s',sVector(2));
                    obj.VectorTex   = sprintf('E_{%s}',sVector(2));
                case 'b'
                    obj.VectorName  = sprintf('%s B-Field',sAxisName);
                    obj.VectorShort = sprintf('B%s',sVector(2));
                    obj.VectorTex   = sprintf('B_{%s}',sVector(2));
                case 'w'
                    obj.VectorName  = sprintf('%s Wakefield',sAxisName);
                    obj.VectorShort = sprintf('W%s',sVector(2));
                    obj.VectorTex   = sprintf('W_{%s}',sVector(2));
                case 'j'
                    obj.VectorName  = sprintf('%s Plasma Current Density',sAxisName);
                    obj.VectorShort = sprintf('Jp%s',sVector(3));
                    obj.VectorTex   = sprintf('J_{p,%s}',sVector(3));
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
                    switch obj.VectorVar
                        case 'wz'
                            aE    = obj.Data.Data(obj.Time,'F','EZ','',sSlice);
                            aData = aE;
                        case 'wr'
                            switch sSlice
                                case 'XZ'
                                    aE    = obj.Data.Data(obj.Time,'F','EX','','XZ');
                                    aB1   = obj.Data.Data(obj.Time,'F','BY','','XZ');
                                    aData = aE - aB1;
                                case 'YZ'
                                    aE    = obj.Data.Data(obj.Time,'F','EY','','YZ');
                                    aB1   = obj.Data.Data(obj.Time,'F','BX','','YZ');
                                    aData = aE + aB1;
                            end % switch
                    end % switch
            end % switch
            if isempty(aData)
                fprintf(2,'QPICVector.Density2D: No data\n');
                return;
            end % if
            
            % Slice and crop
            stData = obj.fParseGridData2D(aData);
            aData  = stData.Data;
            aHAxis = stData.HAxis;
            aVAxis = stData.VAxis;
            cAxes  = stData.Axes;
            
            % Return Data
            stReturn.Data  = aData*obj.VectorFac;
            stReturn.HAxis = aHAxis;
            stReturn.VAxis = aVAxis;
            stReturn.Axes  = {QPICTools.fLabelAxis(cAxes{1}) QPICTools.fLabelAxis(cAxes{2})};
            stReturn.ZPos  = obj.fGetZPos();        
        
        end % function

        function stReturn = Lineout(obj, sSlice, sAxis, aStart, aAverage)

            % Input/Output
            stReturn = {};
            
            if nargin < 5
                aAverage = 1;
            end % if
            
            if nargin < 4
                aStart = 1;
            end % if

            if nargin < 3
                sAxis = 'Z';
            end % if

            if nargin < 2
                sSlice = '';
            end % if
            
            [sAxis,iAxis] = QPICTools.fTranslateAxis(sAxis);
            
            % Get Data and Parse it
            switch obj.VectorType
                case 'e'
                    aData = obj.Data.Data(obj.Time,'F',obj.VectorVar,'',sSlice);
                case 'b'
                    aData = obj.Data.Data(obj.Time,'F',obj.VectorVar,'',sSlice);
                case 'j'
                    aData = obj.Data.Data(obj.Time,'J',obj.VectorVar,'',sSlice);
                case 'w'
                    fprintf(2,'QPICVector.Lineout: For wakefields use WFLineout instead.\n');
                    return;
            end % switch
            if isempty(aData)
                fprintf(2,'QPICVector.Lineout: No data\n');
                return;
            end % if

            stData = obj.fParseGridData1D(aData,sSlice,iAxis,aStart,aAverage);
            if isempty(stData)
                return;
            end % if

            % Return Data
            stReturn.Data   = stData.Data*obj.ScalarFac;
            stReturn.HAxis  = stData.HAxis;
            stReturn.HRange = stData.HLim;
            stReturn.Axis   = sAxis;
            stReturn.ZPos   = obj.fGetZPos();        
        
        end % function
        
        function stReturn = Wakefield2D(obj, sSlice)
            
            % Input/Output
            stReturn = {};
                        
        end % function

        function stReturn = WFLineout(obj, sSlice, sAxis, aStart, aAverage)

            % Input/Output
            stReturn = {};
            
            if nargin < 5
                aAverage = 1;
            end % if
            
            if nargin < 4
                aStart = 1;
            end % if

            if nargin < 3
                sAxis = 'Z';
            end % if

            if nargin < 2
                sSlice = '';
            end % if
            
            [sAxis,iAxis] = QPICTools.fTranslateAxis(sAxis);
            
            if obj.VectorType ~= 'w'
                fprintf(2,'QPICVector.WFLineout: For non-wakefields use Lineout instead.\n');
                return;
            end % if
            
            switch obj.VectorVar
            end % switch

            % Get Data and Parse it
            stData = obj.fParseGridData1D(aData,sSlice,iAxis,aStart,aAverage);
            if isempty(stData)
                return;
            end % if

            % Return Data
            stReturn.Data   = stData.Data*obj.ScalarFac;
            stReturn.HAxis  = stData.HAxis;
            stReturn.HRange = stData.HLim;
            stReturn.Axis   = sAxis;
            stReturn.ZPos   = obj.fGetZPos();        
          
        end % function
        
        
    end % methods

end % classdef
