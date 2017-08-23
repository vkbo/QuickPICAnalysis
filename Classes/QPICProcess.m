
%
%  QPICProcess
% *************
%  A set of static functions for QuickPIC Analysis
%  The classdef only serves as a namespace
%

classdef QPICProcess

    methods(Static, Access='public')
        
        function stReturn = fDeposit(aData, aWeights, iOrder, aGrid, aLim)
            
            % Input/Output
            stReturn = {};
            
            if nargin == 0
                fprintf('\n');
                fprintf('  Function: fDeposit\n');
                fprintf(' ********************\n');
                fprintf('  Deposits weighted data onto grid.\n');
                fprintf('\n');
                fprintf('  Inputs:\n');
                fprintf(' =========\n');
                fprintf('  aData    :: Data coordinates, one column for each dimension.\n');
                fprintf('  aWeights :: One column of weights.\n');
                fprintf('  iOrder   :: Interpolation order. 0 or 1.\n');
                fprintf('  aGrid    :: Vector of resolution (bins), one for each dimension.\n');
                fprintf('  aLim     :: Limits. If left out, uses min and max of each data vector.\n');
                fprintf('\n');
                return;
            end % if

            if nargin < 5
                aLim = [];
            end % if
            
            if iOrder < 0 || iOrder > 1
                fprintf(2,'Error: Order must be 0 or 1.\n');
                return;
            end % if
            
            [iN, iDim]   = size(aData);
            [iNW, iDimW] = size(aWeights);
            iDimG        = numel(aGrid);
            iDimL        = numel(aLim);
            aLim         = aLim(:);
            
            if iDimW ~= 1 || iN ~= iNW || iDimG ~= iDim || iDimL ~= 0 && iDimL ~= 2*iDim
                fprintf(2,'Error: Dimension mismatch in input data.\n');
                return;
            end % if
            
            if isempty(aLim)
                aLim = zeros(2*iDim,1);
                for d=1:iDim
                    aLim(2*d-1) = min(aData(:,d));
                    aLim(2*d)   = max(aData(:,d));
                end % for
                aLim = reshape(aLim,2,iDim);
            else % Prune dataset
                aLim = reshape(aLim,2,iDim);
                for d=1:iDim
                    aData(aData(:,d) < aLim(1,d) | aData(:,d) > aLim(2,d),:) = [];
                end % for
                [iN,~] = size(aData);
            end % if
            aSpan = aLim(2,:)-aLim(1,:);
            aData = aData-aLim(1,:);
            
            aDX = zeros(1,iDim);
            for d=1:iDim
                aDX(d) = aSpan(d)/(aGrid(d)-1);
            end % for
            
            %
            % These are separated by dimension and order to optimise speed for each
            %
            
            % 1D, 0th Order
            if iDim == 1 && iOrder == 0

                aDeposit = zeros(aGrid,1);
                dDX = aDX(1);

                for n=1:iN
                    iPos = round(aData(n)/dDX)+1;
                    aDeposit(iPos) = aDeposit(iPos) + aWeights(n);
                end % for

            end % if
            
            % 2D, 0th Order
            if iDim == 2 && iOrder == 0

                aDeposit = zeros(aGrid);
                dDX1 = aDX(1);
                dDX2 = aDX(2);

                for n=1:iN
                    iPos1 = round(aData(n,1)/dDX1)+1;
                    iPos2 = round(aData(n,2)/dDX2)+1;
                    aDeposit(iPos1,iPos2) = aDeposit(iPos1,iPos2) + aWeights(n);
                end % for

            end % if

            % 1D, 1st Order
            if iDim == 1 && iOrder == 1

                aDeposit = zeros(aGrid+1,1);
                dDX = aDX(1);

                for n=1:iN
                    dPos = aData(n)/dDX+1;
                    iPos = floor(dPos);
                    dRem = dPos - iPos;
                    aDeposit(iPos)   = aDeposit(iPos)   + (1-dRem)*aWeights(n);
                    aDeposit(iPos+1) = aDeposit(iPos+1) + dRem*aWeights(n);
                end % for
                aDeposit = aDeposit(1:end-1);

            end % if
           
            % 2D, 1st Order
            if iDim == 2 && iOrder == 1

                aDeposit = zeros(aGrid(1)+1,aGrid(2)+1);
                dDX1 = aDX(1);
                dDX2 = aDX(2);

                for n=1:iN
                    dPos1 = aData(n,1)/dDX1+1;
                    dPos2 = aData(n,2)/dDX2+1;
                    iPos1 = floor(dPos1);
                    iPos2 = floor(dPos2);
                    dRem1 = dPos1 - iPos1;
                    dRem2 = dPos2 - iPos2;
                    aDeposit(iPos1,   iPos2)   = aDeposit(iPos1,   iPos2)   + (1-dRem1) * (1-dRem2) * aWeights(n);
                    aDeposit(iPos1+1, iPos2)   = aDeposit(iPos1+1, iPos2)   +    dRem1  * (1-dRem2) * aWeights(n);
                    aDeposit(iPos1,   iPos2+1) = aDeposit(iPos1,   iPos2+1) + (1-dRem1) *    dRem2  * aWeights(n);
                    aDeposit(iPos1+1, iPos2+1) = aDeposit(iPos1+1, iPos2+1) +    dRem1  *    dRem2  * aWeights(n);
                end % for
                aDeposit = aDeposit(1:end-1,1:end-1);

            end % if
           
            % Return Results
            stReturn.Deposit = aDeposit;
            stReturn.DeltaX  = aDX;
            stReturn.Limits  = aLim;
            
        end % function
        
    end % methods
    
end % classdef
