
%
%  Function: plotQPSliceEmittance
% ********************************
%  Plots emiitance slice by slice for a beam
%
%  Inputs:
% =========
%  oData   :: QPICData object
%  iTime   :: Time dump
%  sBeam   :: Which beam to look at
%  sAxis   :: Which axis, X or Y
%
%  Options:
% ==========
%  Scan        :: Number of ±sigma beam to scan over. Default 5.0
%  Smooth      :: Number of dZ cells to smooth over for emittance calculation
%  MinStat     :: Minimum number of macroparticles to accept for emittance calculation
%  FigureSize  :: Default [900 500]
%  IsSubplot   :: Default No
%  AutoResize  :: Default On
%

function stReturn = plotQPSliceEmittance(oData, iTime, iBeam, sAxis, varargin)

    % Input/Output

    stReturn = {};

    if nargin == 0
        fprintf('\n');
        fprintf('\n');
        return;
    end % if

    oOpt = inputParser;
    addParameter(oOpt, 'Scan',        5.0);
    addParameter(oOpt, 'EmitTol',     5.0);
    addParameter(oOpt, 'MomTol',      3.0);
    addParameter(oOpt, 'Smooth',      4.0);
    addParameter(oOpt, 'MinStat',     100);
    addParameter(oOpt, 'FigureSize',  [1600 600]);
    addParameter(oOpt, 'IsSubPlot',   'No');
    addParameter(oOpt, 'AutoResize',  'On');
    addParameter(oOpt, 'Absolute',    'Yes');
    parse(oOpt, varargin{:});
    stOpt = oOpt.Results;
    
    switch lower(sAxis)
        case 'x'
            sAxis  = 'x';
            sSlice = 'XZ';
        case 'y'
            sAxis  = 'y';
            sSlice = 'YZ';
        otherwise
            fprintf(2,'Error: Invalid axis. Must be x or y.\n');
            return;
    end % switch

    sBeam  = sprintf('Beam%02d',iBeam);

    dLFac  = oData.Config.Convert.SI.LengthFac;
    dDT    = oData.Config.Simulation.TimeStep;
    dZMax  = oData.Config.Simulation.XMax(1);
    aGrid  = oData.Config.Simulation.Grid;
    dPos   = oData.Config.Beam.(sBeam).Position(1)/dLFac;
    dSz    = oData.Config.Beam.(sBeam).Profile(1)/dLFac;
    dDz    = dZMax/double(aGrid(1))/dLFac;
    
    % Get Beam and Plasma Lineouts

    aStart = aGrid(2:3)/2-1;
    aAvg   = [2 2];
    
    oEB    = QPICScalar(oData,'EB','Units','SI','Scale','mm');
    oEB.Time = iTime;

    stData = oEB.Lineout(sSlice,'Z',aStart,aAvg);
    aQEB   = stData.Data;
    aAxis  = stData.HAxis;

    oEP    = QPICScalar(oData,'EP01','Units','SI','Scale','mm');
    oEP.Time = iTime;

    stData = oEP.Lineout(sSlice,'Z',aStart,aAvg);
    aQEP   = stData.Data;
    
    if strcmpi(stOpt.Absolute, 'Yes')
        aQEB = abs(aQEB);
        aQEP = abs(aQEP);
    end % if
    
    % Crop Datasets
    iMinZ  = round((dPos - dSz*stOpt.Scan)/dDz);
    iMaxZ  = round((dPos + dSz*stOpt.Scan)/dDz);

    aQEB   = aQEB(iMinZ:iMaxZ);
    aQEP   = aQEP(iMinZ:iMaxZ);
    aAxis  = aAxis(iMinZ:iMaxZ);

    dPFac  = 10^floor(log10(max(aQEB)));
    aQEP   = aQEP*dPFac;
    
    % Emittance
    
    oBeam  = QPICBeam(oData,iBeam,'Units','SI','Scale','mm');
    oBeam.Time = iTime;

    stData = oBeam.SlicedPhaseSpace('Dimension',sAxis,'Lim',[iMinZ iMaxZ], ...
                                    'Smooth',stOpt.Smooth,'MinStat',stOpt.MinStat, ...
                                    'EmitTol',stOpt.EmitTol,'ReturnInc','Yes');
    
    aEmitt = stData.ENorm;
    aExcl  = stData.Excluded;
    aMom   = stData.Momentum;
    dMMom  = stData.MeanMom;
    dEMax  = stData.ETolerance;
    dQTot  = stData.TotCharge;
    dQInc  = stData.TolCharge;
    
    % Momentum
    
    aInc   = stData.Included;
    stMom  = oBeam.ScanMomentum(aInc,'Resolution',1e6,'Tolerance',stOpt.MomTol);
    
    aPAxis = stMom.PValues; % Axis with momentum values
    aGAxis = stMom.GValues; % Axis with gamma values
    aQTot  = stMom.Charge;  % Charge per interval
    
    
    % Make a historgam of the P distribution with the same limits for reference
    aHLim  = [aGAxis(1) aGAxis(end)];
    iHNum  = numel(aGAxis);
    stHist = QPICProcess.fDeposit(aInc(:,4),ones(numel(aInc(:,4)),1),1,iHNum,aHLim);
    aHist  = stHist.Deposit;
    
    
    % Scale Data
    
    dMMax  = max(abs(aMom));
    [dTemp,sMomUnit] = QPICTools.fAutoScale(dMMax,stData.MomUnit);
    aMom   = aMom*dTemp/dMMax;
    dMMom  = dMMom*dTemp/dMMax;
    
    [dTemp,sQUnit] = QPICTools.fAutoScale(dQTot,stData.ChargeUnit);
    dQInc  = dQInc*dTemp/dQTot;
    dQTot  = dTemp;
    dQRat  = dQInc/dQTot;
    
    iExcl  = sum(aExcl);
    if iExcl > 0
        fprintf(2,'Warning: %d slices were excluded from emittance calculation due to low statistics.\n',iExcl);
    end % if


    % Plot
    
    if strcmpi(stOpt.IsSubPlot, 'No')
        clf;
        if strcmpi(stOpt.AutoResize, 'On')
            QPICTools.fFigureSize(gcf, stOpt.FigureSize);
        end % if
        set(gcf,'Name',sprintf('%s (%s #%d)',oEB.ScalarName,oData.Config.Name,iTime))
    else
        cla;
    end % if
    
    % Density and Emittance Plot
    
    subplot(4,2,[1,3,5]);

    hold on;

    yyaxis left;
    plot(aAxis,aQEB);
    plot(aAxis,aQEP);

    ylabel('|N_{b,p}/N_0|');

    yyaxis right;
    stairs(aAxis,aEmitt);
    plot(aAxis,ones(numel(aEmitt),1)*dEMax);
    
    text(aAxis(5),dEMax,sprintf('Lim = %.2f µm\nQ_b = %.2f %s [%.1f %%]',dEMax,dQInc,sQUnit,100*dQRat));

    ylabel('\epsilon_N [µm]');
    
    hold off;
    
    title(sprintf('Charge Density and %s-Emittance at %s',upper(sAxis),oEB.PlasmaPosition));
    
    cLegend{1} = 'N_b';
    cLegend{2} = 'N_p';
    cLegend{3} = '\epsilon_N(\xi)';
    %cLegend{4} = sprintf('Q_b = %.2f %s',dQInc,sQUnit);

    if dPFac > 1.0
        cLegend{2} = sprintf('%d·N_p',dPFac);
    end % if
    if dPFac < 0.0
        cLegend{2} = sprintf('N_p/%d',1/dPFac);
    end % if
    
    legend(cLegend,'Location','NW');

    xlim([aAxis(1) aAxis(end)]);
    
    % Forward Momentum Plot
    
    subplot(4,2,7);
    
    aCol = get(gca,'colororder');

    hold on;
    
    plot(aAxis,aMom,'Color',aCol(1,:));
    plot(aAxis,ones(numel(aMom),1)*dMMom,'Color',aCol(1,:),'LineStyle','--');
    
    hold off;
    
    xlabel('\xi - \mu_{b,0} [µm]');
    ylabel(sprintf('P_z [%s]',sMomUnit));

    xlim([aAxis(1) aAxis(end)]);
    ylim([0.9*min(aMom(aMom > 0)) 1.1*max(aMom)]);
    
    % Momentum Histogram
    
    subplot(4,2,[2,4,6]);
    
    hold on;

    yyaxis left;
    plot(aPAxis,abs(aQTot)*1e12);
    
    ylabel(sprintf('{\\fontsize{15}\\int}_{\\fontsize{6}-%.1f%%}^{\\fontsize{6}+%.1f%%} Q dP [pC]',stOpt.MomTol,stOpt.MomTol));

    yyaxis right;
    stairs(aPAxis,aHist);
    ylabel('N / \DeltaP');
    
    xlim([aPAxis(1) aPAxis(end)]);
    
    hold off;

end % function
