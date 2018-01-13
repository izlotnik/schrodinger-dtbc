function [ pError, time, pVisual, pCalc, sFileName, U, U1, U2, U3, U4 ] = calcSolutionRichardson(pTask, pInitWave, solution, error, plotting)
nM_R = solution('M')*(solution('r'):-1:1)/solution('r');
f = GetValue(error, 'frequency', 1);
r = length(nM_R);
time = containers.Map({'total','solve','lhs','conv','rhs','others','kernel','matrix','error','plot','file'}, {0, 0, 0, 0, 0, 0, [0 0], 0, 0, 0, 0});

%% plot
tic
if plotting('details')
    plotErrorL2 = []; plotErrorL2_R = [];
    plotErrorC  = []; plotErrorC_R  = [];
    error_r = error;
    plotting_r = plotting;
else
    error_r = containers.Map('KeyType', 'char', 'ValueType', 'any'); error_r('calc') = false;
    plotting_r = containers.Map({'solution', 'potential', 'count', 'solution_norm', 'norm', 'details'}, {false, false, 0, false, false, false});
end
time('plot') = time('plot') + toc;

%% solve, plot
for m=1:r
    solution('M') = nM_R(m);
    [ nError, time_r, ~, pCalc, ~ ] = calcSolutionCN(pTask, pInitWave, solution, error_r, plotting_r);
    for t=time_r.keys()
        time(t{:}) = time(t{:}) + time_r(t{:});
    end
    switch m
        case 1, U1 = pCalc.U;
        case 2, U2 = pCalc.U;
        case 3, U3 = pCalc.U;
        case 4, U4 = pCalc.U;
    end
    tic
    if plotting('details')
        plotErrorL2   = [ nError.L2Abs(1:((r-m+1)*f):end); plotErrorL2   ];
        plotErrorL2_R = [ nError.L2Rel(1:((r-m+1)*f):end); plotErrorL2_R ];
        plotErrorC    = [ nError.CAbs( 1:((r-m+1)*f):end); plotErrorC    ];
        plotErrorC_R  = [ nError.CRel( 1:((r-m+1)*f):end); plotErrorC_R  ];
    end
    time('plot') = time('plot') + toc;
end

%% solve
tic
nK_R = {1, [ 4/3 -1/3 ], [ 81/40 -16/15 1/24 ], [ 1024/315 -729/280 16/45 -1/360 ]};
switch r
    case 1
        U = nK_R{r}(1) * U1(:, 1:r:end);
        U2 = []; U3 = []; U4 = [];
    case 2
        U = nK_R{r}(1) * U1(:, 1:r:end) + nK_R{r}(2) * U2(:, 1:(r-1):end);
        U3 = []; U4 = [];
    case 3
        U = nK_R{r}(1) * U1(:, 1:r:end) + nK_R{r}(2) * U2(:, 1:(r-1):end) + nK_R{r}(3) * U3(:, 1:(r-2):end);
        U4 = [];
    case 4
        U = nK_R{r}(1) * U1(:, 1:r:end) + nK_R{r}(2) * U2(:, 1:(r-1):end) + nK_R{r}(3) * U3(:, 1:(r-2):end) + nK_R{r}(4) * U4(:, 1:(r-3):end);
end
time('others') = toc;

pCalc.U = U;

%% plot
tic
% if plotting('solution') || plotting('norm')
[pVisual, pCalc] = setVisualParams(pTask, pCalc, ...
    [ plotting('solution') plotting('potential') plotting('count') 0 0 ], ...
    [ 1 f 1 0 ], [0 1]);
% if isfield(pVisual.pPlot2D, 'nLimitY')
%     pVisual.pPlot2D = rmfield(pVisual.pPlot2D, 'nLimitY');
% end
if isfield(pVisual.pPlotError, 'sLabelY')
    pVisual.pPlotError = rmfield(pVisual.pPlotError, 'sLabelY');
end
sFileName = [ pVisual.pPlotError.sPath '/' pVisual.pPlotError.sFilePrefix '_' pVisual.pPlotError.sFileName '_' pVisual.pPlotError.sFileSuffix '_' 'RICH(' num2str(r) ')' ];
% end
if plotting('solution')
    pVisual.pPlot2D = initGraphics('2D', pVisual.pPlot2D);
    for k=1:pCalc.m
        if sum(k-1 == linspace(0, pCalc.m-1, plotting('count')+1))>0
            switch pCalc.name
                case 'FFDS'
                    Ua = U(pCalc.n_perInner, k);
                case 'FEM'
                    Ua = U(pCalc.n_per     , k);
            end
            pVisual.pPlot2D = plot2DSolutionBehaviour(Ua, pVisual.pPlot2D, pTask, k-1);
            saveImage(pVisual.pPlot2D.hFigure, pVisual.pPlot2D, [ 'k=' num2str(k-1) ], mfilename);
        end
    end
    close(gcf);
end
time('plot') = time('plot') + toc;

%% error
tic
if error('calc')
    pError = calcErrorsAll(U, pTask, pInitWave, pCalc, error);
    pError.CAbs(1) = 0; pError.L2Abs(1) = 0;
    pError.CRel(1) = 0; pError.L2Rel(1) = 0;
    printMaxErrors(pError, [], pCalc);
    if plotting('norm')
        pVisual.pPlotError.sTitle = {'', ''};
        hf = plotErrorBehaviour(pVisual.pPlotError, ...
            [ pError.CAbs(pCalc.nFreq); pError.L2Abs(pCalc.nFreq) ], ...
            [ pError.CRel(pCalc.nFreq); pError.L2Rel(pCalc.nFreq) ]);
        saveImage(hf, pVisual.pPlotError, [ 'RICH(' num2str(r) ')' ], mfilename); close(hf);
    end
    if plotting('details')
        p = plotInit;
        switch r
            case 1
                sLineStyle   = {'--', p.line{mod(r+10-1, length(p.line))+1}};
                sColor       = {'b', p.color(mod(r+10-1, length(p.color ))+1, :)};
                sLegendTitle = {'$k=1$', '$r=1$'};
            case 2
                sLineStyle   = {'-', '--', p.line{mod(r+10-1, length(p.line))+1}};
                sColor       = {'r', 'b', p.color(mod(r+10-1, length(p.color))+1, :)};
                sLegendTitle = {'$k=2$', '$k=1$', '$r=2$'};
            case 3
                sLineStyle   = {'-', '-.', '--', p.line{mod(r+10-1, length(p.line))+1}};
                sColor       = {'k', 'g', 'b', p.color(mod(r+10-1, length(p.color))+1, :)};
                sLegendTitle = {'$k=3$', '$k=3/2$', '$k=1$', '$r=3$'};
            case 4
                sLineStyle   = {'--', '-', '-.', '--', p.line{mod(r+10-1, length(p.line))+1}};
                sColor       = {'m', 'y', 'c', 'b', p.color(mod(r+10-1, length(p.color))+1, :)};
                sLegendTitle = {'$k=4$', '$k=2$', '$k=4/3$', '$k=1$', '$r=4$'};
        end
        pVisual.pPlotError.sLineStyle       = sLineStyle;
        pVisual.pPlotError.sColor           = sColor;
        pVisual.pPlotError.sYScale          = 'log'; % 'linear'
        pVisual.pPlotError.sLegendLocation  = 'West'; % 'East'; % 'Best';
        pVisual.pPlotError.sLegendTitle     = sLegendTitle;
        pVisual.pPlotError.sTitle           = {'', ''}; % {'L2-abs', 'C-abs'};
        pVisual.pPlotError.bAlignYLim       = false; % true; %
        pVisual.pPlotError.nLimY            = [10^(-10) 10^(-1)]; % EX22m
        hf = plotErrorBehaviour(pVisual.pPlotError, [ plotErrorL2; pError.L2Abs(1:f:end) ], [ plotErrorC; pError.CAbs(1:f:end) ]);
        saveImage(hf, pVisual.pPlotError, [ 'RICH(' num2str(r) ')_DETAILS' ], mfilename); close(hf);
        if false % Plot relative errors
            pVisual.pPlotError.sTitle       = {'', ''}; % {'L2-rel', 'C-rel'};
            hf = plotErrorBehaviour(pVisual.pPlotError, [ plotErrorL2_R; pError.L2Rel(1:nErrorFrequency:end) ], [ plotErrorC_R; pError.CRel(1:nErrorFrequency:end) ]);
            saveImage(hf, pVisual.pPlotError, [ 'RICH(' num2str(r) ')_DETAILS_REL' ], mfilename); close(hf);
        end
    end
else
    pError = [];
end
time('error') = time('error') + toc;

%%
% time('solve') = time('solve') + (time('lhs')+time('conv')+time('rhs')); time('others') = 0;
% time('total') = time('solve') + sum(time('kernel')) + time('matrix');
% time('others') = time('solve') - (time('lhs')+time('conv')+time('rhs'));
time('solve') = time('solve') + time('others');
time('total') = time('total') + time('others');
T = datestr(datenum(0,0,0,0,0,[ time('total') time('solve') time('lhs') time('conv') time('rhs') ...
    time('others') time('kernel') time('matrix') time('error') time('plot') time('file') ]), 'HH:MM:SS.FFF');
fprintf('%s Timing:\n total=%s\n\t solve =%s\n\t\t lhs   =%s\n\t\t conv  =%s\n\t\t rhs   =%s\n\t\t others=%s\n\t kernel=(%s,%s)\n\t matrix=%s\n error=%s\n plot =%s\n file =%s\n', ...
    datestr(now, 'HH:MM:SS'), T(1, :), T(2, :), T(3, :), T(4, :), T(5, :), T(6, :), T(7, :), T(8, :), T(9, :), T(10, :), T(11, :), T(12, :));