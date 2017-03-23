
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
    addParameter(oOpt, 'Smooth',      4.0);
    addParameter(oOpt, 'MinStat',     100);
    addParameter(oOpt, 'FigureSize',  [900 500]);
    addParameter(oOpt, 'IsSubPlot',   'No');
    addParameter(oOpt, 'AutoResize',  'On');
    addParameter(oOpt, 'Absolute',    'Yes');
    parse(oOpt, varargin{:});
    stOpt = oOpt.Results;
    
    switch lower(sAxis)
        case 'x'
            sAxis  = 'x';
            sSlice = 'XZ';
            iXDim  = 2;
            iPDim  = 5;
        case 'y'
            sAxis  = 'y';
            sSlice = 'YZ';
            iXDim  = 3;
            iPDim  = 6;
        otherwise
            fprintf(2,'Error: Invalid axis. Must be x or y.\n');
            return;
    end % switch

    sBeam  = sprintf('Beam%02d',iBeam);
    dScan  = stOpt.Scan;
    dAvg   = stOpt.Smooth;
    iMinS  = stOpt.MinStat;

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
    iMinZ  = round((dPos - dSz*dScan)/dDz);
    iMaxZ  = round((dPos + dSz*dScan)/dDz);

    aQEB   = aQEB(iMinZ:iMaxZ);
    aQEP   = aQEP(iMinZ:iMaxZ);
    aAxis  = aAxis(iMinZ:iMaxZ);

    dPFac  = 10^floor(log10(max(aQEB)));
    aQEP   = aQEP*dPFac;
    
    % Emittance
    
    oBeam  = QPICBeam(oData,iBeam,'Units','SI','Scale','mm');
    stData = oBeam.SlicedPhaseSpace('Dimension',sAxis,'Lim',[iMinZ iMaxZ],'Smooth',stOpt.Smooth,'MinStat',stOpt.MinStat);
    
    aEmitt = stData.ENorm;
    aExcl  = stData.Excluded;
    
    iExcl = sum(aExcl);
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


    hold on;

    yyaxis left;
    plot(aAxis,aQEB);
    plot(aAxis,aQEP);

    xlabel('\xi - \mu_{b,0} [µm]');
    ylabel(sprintf('|N_b/N_0| & |%d \\times N_p/N_0|',dPFac));

    yyaxis right;
    stairs(aAxis,aEmitt);
    ylabel('\epsilon_N [µm]');
    
    hold off;
    
    title(sprintf('Charge Density and %s-Emittance at %s',upper(sAxis),oEB.PlasmaPosition));
    
    cLegend{1} = 'N_b';
    cLegend{2} = 'N_p';        %sprintf('N_p \\times %d',dPFac);
    cLegend{3} = '\epsilon_N'; %sprintf('\\epsilon_N(\\xi \\pm %.1f\\times\\deltaz)',0.5*dAvg);
    legend(cLegend,'Location','NW');

    xlim([aAxis(1) aAxis(end)]);

end % function
