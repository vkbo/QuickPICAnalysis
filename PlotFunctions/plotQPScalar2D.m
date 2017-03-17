
%
%  Function: plotQPScalar2D
% **************************
%  Plots scalar data from QuickPIC in 2D
%
%  Inputs:
% =========
%  oData   :: QPICData object
%  iTime   :: Time dump
%  sScalar :: Which scalar data to look at
%  sSlice  :: Which slice to look at
%
%  Options:
% ==========
%  Limits      :: Axis limits
%  Slice       :: 2D slice coordinate for 3D data
%  SliceAxis   :: 2D slice axis for 3D data
%  FigureSize  :: Default [900 500]
%  HideDump    :: Default No
%  IsSubplot   :: Default No
%  AutoResize  :: Default On
%  CAxis       :: Color axis limits
%  ShowOverlay :: Default Yes
%  Absolute    :: Use absolute values. Default No
%

function stReturn = plotQPScalar2D(oData, iTime, sScalar, sSlice, varargin)

    % Input/Output

    stReturn = {};

    if nargin == 0
        fprintf('\n');
        fprintf('  Function: plotQPScalar2D\n');
        fprintf(' **************************\n');
        fprintf('  Plots scalar data from QuickPIC in 2D\n');
        fprintf('\n');
        fprintf('  Inputs:\n');
        fprintf(' =========\n');
        fprintf('  oData   :: QPICData object\n');
        fprintf('  iTime   :: Time dump\n');
        fprintf('  sScalar :: Which scalar data to look at\n');
        fprintf('  sSlice  :: Which slice to look at\n');
        fprintf('\n');
        fprintf('  Options:\n');
        fprintf(' ==========\n');
        fprintf('  Limits      :: Axis limits\n');
        fprintf('  Slice       :: 2D slice coordinate for 3D data\n');
        fprintf('  SliceAxis   :: 2D slice axis for 3D data\n');
        fprintf('  GridDiag    :: Options for grid diagnostics data.\n');
        fprintf('  FigureSize  :: Default [900 500]\n');
        fprintf('  HideDump    :: Default No\n');
        fprintf('  IsSubplot   :: Default No\n');
        fprintf('  AutoResize  :: Default On\n');
        fprintf('  CAxis       :: Color axis limits\n');
        fprintf('  ShowOverlay :: Default Yes\n');
        fprintf('  Absolute    :: Use absolute values. Default No\n');
        fprintf('\n');
        return;
    end % if

    oOpt = inputParser;
    addParameter(oOpt, 'Limits',      []);
    addParameter(oOpt, 'Slice',       0.0);
    addParameter(oOpt, 'SliceAxis',   3);
    addParameter(oOpt, 'GridDiag',    {});
    addParameter(oOpt, 'FigureSize',  [900 500]);
    addParameter(oOpt, 'HideDump',    'No');
    addParameter(oOpt, 'IsSubPlot',   'No');
    addParameter(oOpt, 'AutoResize',  'On');
    addParameter(oOpt, 'CAxis',       []);
    addParameter(oOpt, 'ShowOverlay', 'Yes');
    addParameter(oOpt, 'Absolute',    'No');
    parse(oOpt, varargin{:});
    stOpt = oOpt.Results;

    if ~isempty(stOpt.Limits) && length(stOpt.Limits) ~= 4
        fprintf(2, 'Error: Limits specified, but must be of dimension 4.\n');
        return;
    end % if
    
    % Prepare Data

    oSca           = QPICScalar(oData,sScalar,'Units','SI','Scale','mm','Symmetric','Yes');
    oSca.Time      = iTime;
    oSca.SliceAxis = stOpt.SliceAxis;
    oSca.Slice     = stOpt.Slice;
    
    if length(stOpt.Limits) == 4
        oSca.X1Lim = stOpt.Limits(1:2);
        oSca.X2Lim = stOpt.Limits(3:4);
    end % if
    
    stData = oSca.Density2D(sSlice);

    if isempty(stData)
        fprintf(2, 'Error: No data.\n');
        stReturn.Error = 'No data';
        return;
    end % if

    aData  = stData.Data;
    aHAxis = stData.HAxis;
    aVAxis = stData.VAxis;
    sHAxis = stData.Axes{1};
    sVAxis = stData.Axes{2};
    dZPos  = stData.ZPos;

    dPeak  = max(abs(aData(:)));
    [dTemp, sScaUnit] = QPICTools.fAutoScale(dPeak, oSca.ScalarUnit);
    dScale = dTemp/dPeak;

    stReturn.HAxis     = stData.HAxis;
    stReturn.VAxis     = stData.VAxis;
    stReturn.ZPos      = stData.ZPos;
    stReturn.AxisFac   = oSca.AxisFac;
    stReturn.AxisRange = oSca.AxisRange;
    
    if strcmpi(stOpt.Absolute, 'Yes')
        aData = abs(aData);
    end % if

    % Plot
    
    if strcmpi(stOpt.IsSubPlot, 'No')
        clf;
        if strcmpi(stOpt.AutoResize, 'On')
            QPICTools.fFigureSize(gcf, stOpt.FigureSize);
        end % if
        set(gcf,'Name',sprintf('%s (%s #%d)',oSca.ScalarName,oData.Config.Name,iTime))
    else
        cla;
    end % if

    imagesc(aHAxis, aVAxis, aData*dScale);
    set(gca,'YDir','Normal');
    colormap('hot');
    hCol = colorbar();
    if ~isempty(stOpt.CAxis)
        caxis(stOpt.CAxis);
    end % if

    hold on;

    title(sprintf('%s at %s, %s',oSca.ScalarName,oSca.PlasmaPosition,oSca.SlicePosition(sSlice)));
    xlabel(sprintf('%s [mm]',sHAxis));
    ylabel(sprintf('%s [mm]',sVAxis));
    ylabel(hCol,sprintf('%s [%s]',oSca.ScalarTex,sScaUnit));
    
    hold off;
    
    
    % Return

    stReturn.Field = oSca.ScalarName;
    stReturn.XLim  = xlim;
    stReturn.YLim  = ylim;
    stReturn.CLim  = caxis;

end % function
