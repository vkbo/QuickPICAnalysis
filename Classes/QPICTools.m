
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
                [~, sUnit] = fAutoScale(dMin,sBaseUnit);
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
        
        function sReturn = fCheckSlice(sSlice)
            
            switch(lower(sSlice))
                case 'xy'
                    sReturn = 'XY';
                case 'yx'
                    sReturn = 'XY';
                case 'xz'
                    sReturn = 'XZ';
                case 'zx'
                    sReturn = 'XZ';
                case 'yz'
                    sReturn = 'YZ';
                case 'zy'
                    sReturn = 'YZ';
                otherwise
                    sReturn = '';
            end % switch
            
        end % function

    end % methods
    
end % classdef

