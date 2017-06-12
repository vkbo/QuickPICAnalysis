
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
    addParameter(oOpt, 'Scan',        5.0);        % Sigma beam to scan over
    addParameter(oOpt, 'EmitTol',     5.0);        % Emittance tolerance for inclusion (slice width)
    addParameter(oOpt, 'Smooth',      4.0);        % Delta_z to smooth over for better emittance stats
    addParameter(oOpt, 'MinStat',     100);        % Minimum macroparticles for emittance in slice calculation
    addParameter(oOpt, 'MomTol',      3.0);        % Momentum tolerance for charge calculation (%)
    addParameter(oOpt, 'MomRes',      1e6);        % Momentum scan resolution (MeV)
    addParameter(oOpt, 'PruneTail',   100);        % Particles to remove from momentum tail
    addParameter(oOpt, 'FigureSize',  [1200 600]);
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
    
    %
    % Get Beam and Plasma Lineouts
    %

    aStart = aGrid(2:3)/2-1;
    aAvg   = [2 2];
    
    oEB    = QPICScalar(oData,'EB','Units','SI','Scale','mm');
    oEB.Time = iTime;

    stData = oEB.Lineout(sSlice,'Z',aStart,aAvg);
    aQEB   = stData.Data;
    aZAxis = stData.HAxis;

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
    aZAxis = aZAxis(iMinZ:iMaxZ);

    dNFac  = 10^floor(log10(max(aQEB)));
    aQEP   = aQEP*dNFac;
    
    %
    % Emittance
    %
    
    oBeam  = QPICBeam(oData,iBeam,'Units','SI','Scale','mm');
    oBeam.Time = iTime;

    stData = oBeam.SlicedPhaseSpace('Dimension',sAxis,'Lim',[iMinZ iMaxZ], ...
                                    'Smooth',stOpt.Smooth,'MinStat',stOpt.MinStat, ...
                                    'EmitTol',stOpt.EmitTol,'ReturnInc','Yes');
    
    aEmitt = stData.ENorm;
    aExcl  = stData.Excluded;
    aMom   = stData.Momentum;
    dMMom  = stData.MeanMom;
    dSMom  = stData.IncRMS;
    dEMax  = stData.ETolerance;
    dQTot  = stData.TotCharge;
    dQInc  = stData.TolCharge;
    dDz    = stData.Resolution;
    
    % Scale Data
    dMMax  = max(abs(aMom));
    [dTemp,sMomUnit] = QPICTools.fAutoScale(dMMax,stData.MomUnit);
    aMom   = aMom*dTemp/dMMax;
    dMMom  = dMMom*dTemp/dMMax;
    [dSMom,sSMomUnit] = QPICTools.fAutoScale(dSMom,stData.MomUnit);
    
    [dTemp,sQUnit] = QPICTools.fAutoScale(dQTot,stData.ChargeUnit);
    dQFac  = dTemp/dQTot;
    dQInc  = dQInc*dQFac;
    dQTot  = dTemp;
    dQRat  = dQInc/dQTot;

    [dDz,sDUnit] = QPICTools.fAutoScale(dDz,'m');
    
    iExcl  = sum(aExcl);
    if iExcl > 0
        fprintf(2,'Warning: %d slices were excluded from emittance calculation due to low statistics\n',iExcl);
    end % if

    stReturn.EmittLim     = dEMax;
    stReturn.EmittInclQ   = dQInc/dQFac;
    stReturn.EmittInclRat = dQRat;

    %
    % Momentum
    %
    
    stMom  = oBeam.ScanMomentum(stData.Included,'Resolution',stOpt.MomRes, ...
                                'Tolerance',stOpt.MomTol,'PruneTail',stOpt.PruneTail);
    
    aPz    = stMom.Particles;
    aPAxis = stMom.PValues;
    aPLims = stMom.PIntervals;
    aQTot  = stMom.Charge;

    aPLim  = [aPAxis(1) aPAxis(end)];
    iPNum  = numel(aPAxis);

    dPMax  = max(aPAxis);
    [dTemp,sPUnit] = QPICTools.fAutoScale(dPMax,'eV/c');
    dPFac  = dTemp/dPMax;
    aPAxis = aPAxis*dPFac;
    
    dSign  = sum(abs(aQTot))/sum(aQTot);
    [dQMax,iQMax] = max(abs(aQTot));
    dQMax  = dQMax*dSign*dQFac;
    dPMax  = aPAxis(iQMax);
    aPInt  = aPLims(iQMax,:)*dPFac;
    dPRat  = dQMax/dQTot;
    
    [dPRes,sRUnit] = QPICTools.fAutoScale(stOpt.MomRes,'eV/c');
    
    % Make a historgam of the P distribution with the same limits for reference
    stHist = QPICProcess.fDeposit(aPz,ones(numel(aPz),1),1,iPNum,aPLim);
    aHist  = stHist.Deposit;
    aHLim  = stHist.Limits;
    aHAxis = linspace(aHLim(1),aHLim(2),numel(aHist))*dPFac;

    stReturn.MomBest    = dPMax/dPFac;
    stReturn.MomInclQ   = dQMax/dQFac;
    stReturn.MomInclRat = dPRat;

    %
    % Plot
    %
    
    fMain = gcf; clf;
    if strcmpi(stOpt.AutoResize, 'On')
        QPICTools.fFigureSize(fMain, stOpt.FigureSize);
    end % if
    fMain.Name  = sprintf('Beam Quality (%s #%d)',oData.Config.Name,iTime);
    fMain.Units = 'Pixels';
    aSize       = [fMain.Position(3:4) fMain.Position(3:4)];
    
    aDN = [0.00 0.30 0.50 0.70].*aSize + [70 40 -140 -75];
    aDP = [0.00 0.00 0.50 0.30].*aSize + [70 60 -140 -70];
    aPZ = [0.50 0.30 0.50 0.70].*aSize + [70 40 -140 -75];
    aIN = [0.50 0.00 0.50 0.30].*aSize + [70 40 -140 -50];
    
    axDN = axes('Units','Pixels','Position',aDN);
    axDP = axes('Units','Pixels','Position',aDP);
    axPZ = axes('Units','Pixels','Position',aPZ);
    
    cEmittLim   = {sprintf('%0.2f',dEMax)     'µm'};
    cChargeIncZ = {sprintf('%0.2f',dQInc)     sQUnit};
    cChargeRatZ = {sprintf('%0.2f',100*dQRat) '%'};
    cScanResZ   = {sprintf('%0.2f',dDz)       sDUnit};

    cBestMom    = {sprintf('%0.2f',dPMax)     sMomUnit};
    cRMSMom     = {sprintf('%0.2f',dSMom)     sSMomUnit};
    cChargeIncP = {sprintf('%0.2f',dQMax)     sQUnit};
    cChargeRatP = {sprintf('%0.2f',100*dPRat) '%'};
    cScanResP   = {sprintf('%0.2f',dPRes)     sRUnit};
    
    bgInfo = uibuttongroup('Title','Details','Units','Pixels','Position',aIN);
    
    iY = aIN(4)-45;
    uicontrol(bgInfo,'Style','Text','String','Emittance Limit:', 'Position',[ 10 iY    120 20],'HorizontalAlignment','Left');
    uicontrol(bgInfo,'Style','Text','String',cEmittLim{1},       'Position',[130 iY     50 20],'HorizontalAlignment','Right');
    uicontrol(bgInfo,'Style','Text','String',cEmittLim{2},       'Position',[180 iY     50 20],'HorizontalAlignment','Left');
    
    uicontrol(bgInfo,'Style','Text','String','Charge Included:', 'Position',[ 10 iY-40 120 20],'HorizontalAlignment','Left');
    uicontrol(bgInfo,'Style','Text','String',cChargeIncZ{1},     'Position',[130 iY-40  50 20],'HorizontalAlignment','Right');
    uicontrol(bgInfo,'Style','Text','String',cChargeIncZ{2},     'Position',[180 iY-40  50 20],'HorizontalAlignment','Left');

    uicontrol(bgInfo,'Style','Text','String','Charge Ratio:',    'Position',[ 10 iY-60 120 20],'HorizontalAlignment','Left');
    uicontrol(bgInfo,'Style','Text','String',cChargeRatZ{1},     'Position',[130 iY-60  50 20],'HorizontalAlignment','Right');
    uicontrol(bgInfo,'Style','Text','String',cChargeRatZ{2},     'Position',[180 iY-60  50 20],'HorizontalAlignment','Left');

    uicontrol(bgInfo,'Style','Text','String','Scan Resolution:', 'Position',[ 10 iY-80 120 20],'HorizontalAlignment','Left');
    uicontrol(bgInfo,'Style','Text','String',cScanResZ{1},       'Position',[130 iY-80  50 20],'HorizontalAlignment','Right');
    uicontrol(bgInfo,'Style','Text','String',cScanResZ{2},       'Position',[180 iY-80  50 20],'HorizontalAlignment','Left');

    uicontrol(bgInfo,'Style','Text','String','Best Momentum:',   'Position',[230 iY    120 20],'HorizontalAlignment','Left');
    uicontrol(bgInfo,'Style','Text','String',cBestMom{1},        'Position',[350 iY     50 20],'HorizontalAlignment','Right');
    uicontrol(bgInfo,'Style','Text','String',cBestMom{2},        'Position',[400 iY     50 20],'HorizontalAlignment','Left');
    
    uicontrol(bgInfo,'Style','Text','String','RMS Momentum:',    'Position',[230 iY-20 120 20],'HorizontalAlignment','Left');
    uicontrol(bgInfo,'Style','Text','String',cRMSMom{1},         'Position',[350 iY-20  50 20],'HorizontalAlignment','Right');
    uicontrol(bgInfo,'Style','Text','String',cRMSMom{2},         'Position',[400 iY-20  50 20],'HorizontalAlignment','Left');
    
    uicontrol(bgInfo,'Style','Text','String','Charge Included:', 'Position',[230 iY-40 120 20],'HorizontalAlignment','Left');
    uicontrol(bgInfo,'Style','Text','String',cChargeIncP{1},     'Position',[350 iY-40  50 20],'HorizontalAlignment','Right');
    uicontrol(bgInfo,'Style','Text','String',cChargeIncP{2},     'Position',[400 iY-40  50 20],'HorizontalAlignment','Left');
    
    uicontrol(bgInfo,'Style','Text','String','Charge Ratio:',    'Position',[230 iY-60 120 20],'HorizontalAlignment','Left');
    uicontrol(bgInfo,'Style','Text','String',cChargeRatP{1},     'Position',[350 iY-60  50 20],'HorizontalAlignment','Right');
    uicontrol(bgInfo,'Style','Text','String',cChargeRatP{2},     'Position',[400 iY-60  50 20],'HorizontalAlignment','Left');
    
    uicontrol(bgInfo,'Style','Text','String','Scan Resolution:', 'Position',[230 iY-80 120 20],'HorizontalAlignment','Left');
    uicontrol(bgInfo,'Style','Text','String',cScanResP{1},       'Position',[350 iY-80  50 20],'HorizontalAlignment','Right');
    uicontrol(bgInfo,'Style','Text','String',cScanResP{2},       'Position',[400 iY-80  50 20],'HorizontalAlignment','Left');

    % Density and Emittance Plot
    
    axes(axDN);

    hold on;

    yyaxis left;
    plot(aZAxis,aQEB);
    plot(aZAxis,aQEP);

    ylabel('|N_{b,p}/N_0|');

    yyaxis right;
    stairs(aZAxis,aEmitt);
    plot(aZAxis,ones(numel(aEmitt),1)*dEMax);

    ylabel('\epsilon_N [µm]');
    
    hold off;
    
    title(sprintf('Charge Density and %s-Emittance at %s',upper(sAxis),oEB.PlasmaPosition));
    
    cLegend{1} = 'N_b';
    cLegend{2} = 'N_p';
    cLegend{3} = '\epsilon_N(\xi)';
    
    if log10(dNFac) > 0.0
        cLegend{2} = sprintf('%d·N_p',dNFac);
    end % if
    if log10(dNFac) < 0.0
        cLegend{2} = sprintf('N_p/%d',1/dNFac);
    end % if
    
    legend(cLegend,'Location','NW');

    xlim([aZAxis(1) aZAxis(end)]);
    
    % Forward Momentum Plot
    
    axes(axDP);
    
    aCol = get(gca,'colororder');

    hold on;
    
    plot(aZAxis,aMom,'Color',aCol(1,:));
    plot(aZAxis,ones(numel(aMom),1)*dMMom,'Color',aCol(1,:),'LineStyle','--');
    
    hold off;
    
    %xlabel('\xi - \mu_{b,0} [µm]');
    xlabel('\xi [µm]');
    ylabel(sprintf('P_z [%s]',sMomUnit));

    xlim([aZAxis(1) aZAxis(end)]);
    ylim([0.9*min(aMom(aMom > 0)) 1.1*max(aMom)]);
    
    % Momentum Histogram
    
    axes(axPZ);
    
    hold on;

    yyaxis left;
    plot(aPAxis,abs(aQTot)*1e12);
    ylabel(sprintf('{\\fontsize{15}\\int}_{\\fontsize{6}-%.1f%%}^{\\fontsize{6}+%.1f%%} Q dP [pC]',stOpt.MomTol,stOpt.MomTol));

    yyaxis right;
    stairs(aHAxis,aHist);
    ylabel('N / \DeltaP');
    
    line([1 1]*dPMax   ,ylim,'Color',aCol(1,:),'LineSTyle','--');
    line([1 1]*aPInt(1),ylim,'Color',aCol(1,:),'LineSTyle',':');
    line([1 1]*aPInt(2),ylim,'Color',aCol(1,:),'LineSTyle',':');
    
    xlabel(sprintf('P_z [%s]',sPUnit));
    if aPLim(1) ~= aPLim(2)
        xlim(aPLim*dPFac);
    end % if
    
    hold off;
    
    title(sprintf('Forward Momentum and Integrated Charge for ±%.1f%%',stOpt.MomTol));

    % Return
    stReturn.ZAxis = aZAxis;
    stReturn.PAxis = aPAxis;

end % function
