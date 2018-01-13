function plotMaxValuesFromTable() % TODO: update for 2D data
clearContext;
%% Set params
file = '2D_EX22a_TRISYS_SSP_20121103_200807.csv'; % '1D_EX02m_TRISYS_CN_20121031_160037.csv';
[ ~, ~, tokenindices, ~, ~, ~, ~ ] = regexp(file, '([\w\d]*)\_([\w\d]*)\_(\w*)\_(\w*)\_(\d*\_\d*)\.(\w*)');
ind = tokenindices{:};
sDimension  = file(ind(1, 1):ind(1, 2));
sExampleName= file(ind(2, 1):ind(2, 2));
sMethod     = file(ind(3, 1):ind(3, 2));
sType       = file(ind(4, 1):ind(4, 2));
sSuffix     = file(ind(5, 1):ind(5, 2));

%% Figure params
p.bSaveFigure   = true;
p.sPrintType    = '-depsc2';
p.sPath         = [ cd '\data\plots' ];
p.sFilePrefix   = [ sDimension '_' sExampleName '_' sMethod];
p.sFileName     = sType;
p.sFileSuffix   = sSuffix;
p.sFileExtension= 'eps';

%% Read params
file = [ cd '\data\plots\TableMaxErrors_20121106\' file ];
q = csvread(file);

%% Plot all permutations
X = perms(1:3); % [1 2 3]; % 
Y = 4;          % 'ERROR_ABS_C'
for x=1:size(X, 1)
    for y=length(Y)
        plot_it(q, [X(x, :) Y(y)], 'loglog', p);
    end
end
% [ hf, sSubSuffix ] = plot_it(q, [2 1 3 5], 'semilogx'); saveImage(hf, p, [ sSubSuffix '_' sSuffix ], mfilename); close(hf);
% [ hf, sSubSuffix ] = plot_it(q, [1 2 3 6], 'semilogy'); saveImage(hf, p, [ sSubSuffix '_' sSuffix ], mfilename); close(hf);
% [ hf, sSubSuffix ] = plot_it(q, [2 1 3 7], 'plot'    ); saveImage(hf, p, [ sSubSuffix '_' sSuffix ], mfilename); close(hf);
% [ hf, sSubSuffix ] = plot_it(q, [1 2 3 8], 'loglog'  ); saveImage(hf, p, [ sSubSuffix '_' sSuffix ], mfilename); close(hf);
end

function plot_it(q, nCol, sType, Params)
J = unique(q(:, nCol(1)));
M = unique(q(:, nCol(2)));
K = unique(q(:, nCol(3)));
switch nCol(1)
    case 1
        sLegend = 'J';
    case 2
        sLegend = 'M';
    case 3
        sLegend = 'K';
    otherwise
        sLegend = '' ; warning('Unknown x data column ID');
end
switch nCol(2)
    case 1
        sLabelX = 'J';
    case 2
        sLabelX = 'M';
    case 3
        sLabelX = 'K';
    otherwise
        sLabelX = '' ; warning('Unknown y(x) data column ID');
end
switch nCol(3)
    case 1
        sFix = 'J';
    case 2
        sFix = 'M';
    case 3
        sFix = 'K';
    otherwise
        sFix = '' ; warning('Unknown fixed data column ID');
end
sLegendLocation = 'Best'; % 'SouthWest';
switch nCol(4)
    case 4
        sLabelY = 'Maximum in time absolute C-error' ; sSuffix = 'ERROR_ABS_C' ;
    case 5
        sLabelY = 'Maximum in time absolute L2-error'; sSuffix = 'ERROR_ABS_L2';
    case 6
        sLabelY = 'Maximum in time relative C-error' ; sSuffix = 'ERROR_REL_C' ;
    case 7
        sLabelY = 'Maximum in time relative L2-error'; sSuffix = 'ERROR_REL_L2';
    case 8
        sLabelY = 'Calculation time, sec.';            sSuffix = 'TIME_SOL'    ;
        sLegendLocation = 'NorthWest';
    otherwise
        sLabelY = ''; sSuffix = ''; warning('Unknown y data column ID');
end
if ~isempty(sLegend)
    sSuffix = [ sSuffix '_' sLegend '(' sLabelX ')' ];
end
sSuffixD = sSuffix; sFileNameD = Params.sFileName;
for k=1:length(K)
    r = q(q(:, nCol(3))==K(k), :);
    sSuffix = [ sSuffixD '_' sFix '=' num2str(K(k)) ];
    hf = figure('Units', 'Normalize', 'OuterPosition', [0 0 1 1]); p = plotInit; sLegendTitle = '';
    for j=1:length(J)
        s = r(r(:, nCol(1))==J(j), :);
        sLegendTitle{j} = [ '$' sLegend '=' num2str(J(j)) '$' ];
        switch sType
            case 'loglog'
                loglog( s(:, nCol(2)), s(:, nCol(4)) ...
                    ,'LineWidth',       4 ...
                    ,'LineStyle',       p.line{  mod(j-1, length(p.line  ))+1}    ...
                    ,'Color',           p.color( mod(j-1, length(p.color ))+1, :) ...
                    ,'Marker',          p.marker{mod(j  , length(p.marker))+1}    ... % first element is 'none' - skip it
                    ,'MarkerFaceColor', p.color( mod(j-1, length(p.color ))+1, :) ...
                    ,'MarkerSize',      8 );
            case 'semilogx'
                semilogx( s(:, nCol(2)), s(:, nCol(4)) ...
                    ,'LineWidth',       4 ...
                    ,'LineStyle',       p.line{  mod(j-1, length(p.line  ))+1}    ...
                    ,'Color',           p.color( mod(j-1, length(p.color ))+1, :) ...
                    ,'Marker',          p.marker{mod(j  , length(p.marker))+1}    ... % first element is 'none' - skip it
                    ,'MarkerFaceColor', p.color( mod(j-1, length(p.color ))+1, :) ...
                    ,'MarkerSize',      8 );
            case 'semilogy'
                semilogy( s(:, nCol(2)), s(:, nCol(4)) ...
                    ,'LineWidth',       4 ...
                    ,'LineStyle',       p.line{  mod(j-1, length(p.line  ))+1}    ...
                    ,'Color',           p.color( mod(j-1, length(p.color ))+1, :) ...
                    ,'Marker',          p.marker{mod(j  , length(p.marker))+1}    ... % first element is 'none' - skip it
                    ,'MarkerFaceColor', p.color( mod(j-1, length(p.color ))+1, :) ...
                    ,'MarkerSize',      8 );
            otherwise
                plot( s(:, nCol(2)), s(:, nCol(4)) ...
                    ,'LineWidth',       4 ...
                    ,'LineStyle',       p.line{  mod(j-1, length(p.line  ))+1}    ...
                    ,'Color',           p.color( mod(j-1, length(p.color ))+1, :) ...
                    ,'Marker',          p.marker{mod(j  , length(p.marker))+1}    ... % first element is 'none' - skip it
                    ,'MarkerFaceColor', p.color( mod(j-1, length(p.color ))+1, :) ...
                    ,'MarkerSize',      8 );
        end
        set(gca, 'NextPlot', 'Add');
    end
    set(gca, 'xlim', [ min(M) max(M) ] );
    set(gca, 'FontSize', 14, 'FontWeight', 'Bold');
    legend(gca, sLegendTitle, 'Location', sLegendLocation, 'Interpreter', 'LaTeX');
    xlabel(sLabelX); ylabel(sLabelY);
    Params.sFileName = [ sFileNameD '_' sSuffix ];
    saveImage(hf, Params, '', mfilename); close(hf);
end
end