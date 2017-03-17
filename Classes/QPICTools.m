
%
%  QPICTools
% ***********
%  A set of static functions for QuickPIC Analysis
%  The classdef only serves as a namespace
%

classdef QPICTools
    
    methods(Static, Access='public')

        %  fAutoScale
        % ************
        %  Autoscaling of units

        function [dValue, sUnit] = fAutoScale(dBaseValue, sBaseUnit, dMin)

            if nargin < 3
                dMin = 0.999e-24;
            end % if

            dValue = dBaseValue;
            sUnit  = sBaseUnit;

            if abs(dBaseValue) < dMin
                [~, sUnit] = QPICTools.fAutoScale(dMin,sBaseUnit);
                dValue = 0.0;
                return;
            end % if

            if abs(dBaseValue) > 1.0

                if abs(dBaseValue) > 1e24
                    dValue = dBaseValue*1e-24;
                    sUnit  = strcat('Y',sBaseUnit);
                    return;
                end % if

                if abs(dBaseValue) > 1e21
                    dValue = dBaseValue*1e-21;
                    sUnit  = strcat('Z',sBaseUnit);
                    return;
                end % if

                if abs(dBaseValue) > 1e18
                    dValue = dBaseValue*1e-18;
                    sUnit  = strcat('E',sBaseUnit);
                    return;
                end % if

                if abs(dBaseValue) > 1e15
                    dValue = dBaseValue*1e-15;
                    sUnit  = strcat('P',sBaseUnit);
                    return;
                end % if

                if abs(dBaseValue) > 1e12
                    dValue = dBaseValue*1e-12;
                    sUnit  = strcat('T',sBaseUnit);
                    return;
                end % if

                if abs(dBaseValue) > 1e9
                    dValue = dBaseValue*1e-9;
                    sUnit  = strcat('G',sBaseUnit);
                    return;
                end % if

                if abs(dBaseValue) > 1e6
                    dValue = dBaseValue*1e-6;
                    sUnit  = strcat('M',sBaseUnit);
                    return;
                end % if

                if abs(dBaseValue) > 1e3
                    dValue = dBaseValue*1e-3;
                    sUnit  = strcat('k',sBaseUnit);
                    return;
                end % if

            else

                if abs(dBaseValue) < 0.999e-21
                    dValue = dBaseValue*1e24;
                    sUnit  = strcat('y',sBaseUnit);
                    return;
                end % if

                if abs(dBaseValue) < 0.999e-18
                    dValue = dBaseValue*1e21;
                    sUnit  = strcat('z',sBaseUnit);
                    return;
                end % if

                if abs(dBaseValue) < 0.999e-15
                    dValue = dBaseValue*1e18;
                    sUnit  = strcat('a',sBaseUnit);
                    return;
                end % if

                if abs(dBaseValue) < 0.999e-12
                    dValue = dBaseValue*1e15;
                    sUnit  = strcat('f',sBaseUnit);
                    return;
                end % if

                if abs(dBaseValue) < 0.999e-9
                    dValue = dBaseValue*1e12;
                    sUnit  = strcat('p',sBaseUnit);
                    return;
                end % if

                if abs(dBaseValue) < 0.999e-6
                    dValue = dBaseValue*1e9;
                    sUnit  = strcat('n',sBaseUnit);
                    return;
                end % if

                if abs(dBaseValue) < 0.999e-3
                    dValue = dBaseValue*1e6;
                    sUnit  = strcat('µ',sBaseUnit);
                    return;
                end % if

                if abs(dBaseValue) < 0.999
                    dValue = dBaseValue*1e3;
                    sUnit  = strcat('m',sBaseUnit);
                    return;
                end % if

            end % if

        end % function

        %  fFigureSize
        % *************
        %  Resizes a figure while preserving position

        function fFigureSize(figIn, aSize)

            if strcmpi(get(figIn, 'WindowStyle'), 'Docked')
                return;
            end % if

            aPosition      = get(figIn, 'Position');
            aPosition(3:4) = aSize;
            set(figIn, 'Position', aPosition);

        end % function
        
        %  fGetIndex
        % ***********
        %  Returns the closest index of a value in a vector

        function iIndex = fGetIndex(aVector, dValue)

            iIndex = 0;

            if isempty(dValue) || isempty(aVector)
                return;
            end % if

            for i=1:length(aVector)-1
                if dValue >= aVector(1) && dValue < aVector(i+1)
                    if dValue-aVector(i) < aVector(i+1)-dValue
                        iIndex = i;
                        return;
                    else
                        iIndex = i+1;
                        return;
                    end % if
                end % if
            end % for

            if dValue >= aVector(end-1)-(aVector(end)-aVector(end-1))/2
                iIndex = length(aVector);
            end % if

            % Matlab matrix index cannot be less than 1
            if iIndex < 1
                iIndex = 1;
            end % if

        end % function
        
        %  fCheckSlice
        % *************
        %  Checks a slice entry and returns a valid QuickPIC Slice diagnostics
        
        function [sReturn,iDim,iOrth] = fCheckSlice(sSlice)
            
            switch(lower(sSlice))
                case 'xy'
                    sReturn = 'XY';
                    iDim    = 2;
                    iOrth   = 1;
                case 'yx'
                    sReturn = 'XY';
                    iDim    = 2;
                    iOrth   = 1;
                case 'xz'
                    sReturn = 'XZ';
                    iDim    = 2;
                    iOrth   = 3;
                case 'zx'
                    sReturn = 'XZ';
                    iDim    = 2;
                    iOrth   = 3;
                case 'yz'
                    sReturn = 'YZ';
                    iDim    = 2;
                    iOrth   = 2;
                case 'zy'
                    sReturn = 'YZ';
                    iDim    = 2;
                    iOrth   = 2;
                otherwise
                    sReturn = '';
                    iDim    = 3;
                    iOrth   = 0;
            end % switch
            
        end % function

        %  fLengthScale
        % **************
        %  Converts from one scale unit to another

        function [dScale, sUnit] = fLengthScale(sToUnit, sFromUnit)

            dScale = 1.0;
            sUnit  = 'm';

            if nargin < 2
                sFromUnit = 'm';
            end % if

            switch lower(sFromUnit)
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

        %  fTranslateAxis
        % ****************
        %  Converts an axis to the internally used axis name and number

        function [sReturn,iReturn] = fTranslateAxis(vAxis)
            
            sAxis = num2str(vAxis);
            
            switch lower(sAxis)
                case '1'
                    sReturn = 'x1';
                    iReturn = 1;
                case '2'
                    sReturn = 'x2';
                    iReturn = 2;
                case '3'
                    sReturn = 'x3';
                    iReturn = 3;
                case '4'
                    sReturn = 'r';
                    iReturn = 4;
                case 'x'
                    sReturn = 'x2';
                    iReturn = 2;
                case 'y'
                    sReturn = 'x3';
                    iReturn = 3;
                case 'z'
                    sReturn = 'x1';
                    iReturn = 1;
                case 'x1'
                    sReturn = 'x1';
                    iReturn = 1;
                case 'x2'
                    sReturn = 'x2';
                    iReturn = 2;
                case 'x3'
                    sReturn = 'x3';
                    iReturn = 3;
                case 'xi'
                    sReturn = 'x1';
                    iReturn = 1;
                case 'r'
                    sReturn = 'r';
                    iReturn = 4;
                otherwise
                    sReturn = sAxis;
                    iReturn = 0;
            end % switch
            
        end % function
        
        %  fLabelAxis
        % ************
        %  Converts an axis name to a label

        function sReturn = fLabelAxis(vAxis, bUseBox)
            
            sAxis = num2str(vAxis);
            
            if nargin < 2
                bUseBox = true;
            end % if

            switch lower(sAxis)
                case '1'
                    if bUseBox
                        sReturn = '\xi';
                    else
                        sReturn = 'z';
                    end % if
                case '2'
                    sReturn = 'x';
                case '3'
                    sReturn = 'y';
                case '4'
                    sReturn = 'r';
                case 'x'
                    sReturn = 'x';
                case 'y'
                    sReturn = 'y';
                case 'z'
                    if bUseBox
                        sReturn = '\xi';
                    else
                        sReturn = 'z';
                    end % if
                case 'r'
                    sReturn = 'r';
                case 'xi'
                    if bUseBox
                        sReturn = '\xi';
                    else
                        sReturn = 'z';
                    end % if
                case 'x1'
                    if bUseBox
                        sReturn = '\xi';
                    else
                        sReturn = 'z';
                    end % if
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

    end % methods
    
end % classdef

