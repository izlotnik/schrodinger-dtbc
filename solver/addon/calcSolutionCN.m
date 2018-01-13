function [ rError, time, pVisual, pCalc, sFileName, U ] = calcSolutionCN(pTask, pInitWave, solution, error, plotting)
%%
switch solution('type')
    case 'FEM'
        param = solution('N'); type = 'n=%d';
    case 'FFDS'
        param = solution('theta'); type = 'theta=%.4g';
end
if solution.isKey('r')
    fprintf(sprintf('%s %s, J=%d, r=%d, M=%d\n', datestr(now, 'HH:MM:SS'), type, solution('J'), solution('r'), solution('M')), param);
else
    fprintf(sprintf('%s %s, J=%d, M=%d\n', datestr(now, 'HH:MM:SS'), type, solution('J'), solution('M')), param);
end

%%
% DTBC: 
[pCalc, pFile] = setCalcParams(pTask, {solution('type')}, [ param solution('J') solution('M') NaN 1 ], {solution('matrix')});
% SD: [pCalc, pFile] = setCalcParams(pTask, {solution('type'), 'FFDS'}, [param solution('J') solution('M') 1/4 1], {solution('matrix')});
[pVisual, pCalc] = setVisualParams(pTask, pCalc, ...
    [ plotting('solution') plotting('potential') plotting('count') 0 0 ], ...
    [ 1 GetValue(error, 'frequency', 1) 1 0 ], [0 1]);
% if isfield(pVisual.pPlot2D, 'nLimitY')
%     pVisual.pPlot2D = rmfield(pVisual.pPlot2D, 'nLimitY');
% end
% EX24: [1 1 50 0 0]
% pVisual.pPlot2D.nLineWidth = 3;
% pVisual.pPlot2D.nLimitY = [-0.55 0.55]; % [-0.21 0.21];
pCalc.bCalcError = error('calc');
if error.isKey('etalon') && ( strcmp(error('etalon'), 'global') || strcmp(error('etalon'), 'local') ) % Calc solution and errors in comparison to the etalon
    MAIN1D
    rError = calcErrorsAll(pCalc.U, pTask, pInitWave, pCalc, error);
    [pVisual.pPlotError.sLabelY, pVisual.pPlotError.sFileSuffix] = prepareParams(pTask, pVisual, error('suffix'));
else
    if exist(pFile.pSolution.sFullName, 'file') % load solution
        fprintf('%s Load solution...\n', datestr(now, 'HH:MM:SS'));
        U = loadVariable(pFile.pSolution.sFullName, true, mfilename); pCalc.U = U; 
        time = containers.Map({'solve', 'kernel', 'matrix', 'rhs', 'lhs', 'conv', 'error', 'plot', 'file'}, {0, [0 0], 0, 0, 0, 0, 0, 0, 0});
        if pCalc.bCalcError && exist(pFile.pError.sFullName, 'file') % load error
            fprintf('%s Load error...\n', datestr(now, 'HH:MM:SS'));
            rError = loadVariable(pFile.pError.sFullName, true, mfilename);
            tmp = max(rError.CAbs);
            if isempty(tmp) || isnan(tmp) || isinf(tmp) || tmp <= 0 % if not OK then calculate solution
                MAIN1D
            end
        else % if not OK then calculate solution
            MAIN1D
        end
    else
        MAIN1D
    end
    % pVisual.pPlotError.sLabelY = ''; pVisual.pPlotError.sFileSuffix = '';
end
close(gcf);
if pCalc.bCalcError
    printMaxErrors(rError, [], pCalc);
end
%
if pVisual.bPlotDTBCKernel && ~isempty(pCalc.pRBoundary.R) % Plot DTBC kernel
    PLOTKernel(pCalc.pRBoundary.R, pCalc, pVisual.pPlotKernel); close(gcf);
end
%
pVisual.bPlotErrorTimeBehaviour = plotting('norm');
% if pVisual.bPlotErrorTimeBehaviour && ~isempty(rError) % Plot and save plot of time-level error behaviour
%     hf = plotErrorBehaviour(pVisual.pPlotError, [ rError.CAbs(pCalc.nFreq); rError.L2Abs(pCalc.nFreq) ]...
%         , [ rError.CRel(pCalc.nFreq); rError.L2Rel(pCalc.nFreq) ]);
%     saveImage(hf, pVisual.pPlotError, '', mfilename); close(hf);
% end
if pVisual.bPlotErrorTimeBehaviour && ~isempty(rError) % Plot and save plot of time-level error behaviour
    cYScale         = {'log'};      % {'linear', 'log'}; %
    cLegendLocation = {'NorthEast'};% {'NorthWest', 'NorthEast', 'SouthWest', 'SouthEast', 'Best', 'BestOutside'}; %
    for i=cYScale
        pVisual.pPlotError.sYScale = i{:};
        for j=cLegendLocation
            pVisual.pPlotError.sLegendLocation = j{:};
            hf = plotErrorBehaviour(pVisual.pPlotError, ...
                [ rError.CAbs(pCalc.nFreq); rError.L2Abs(pCalc.nFreq) ], ...
                [ rError.CRel(pCalc.nFreq); rError.L2Rel(pCalc.nFreq) ]);
            if isfield(pVisual.pPlotError, 'sYScale') && ~strcmp(pVisual.pPlotError.sYScale, 'linear')
                sPrefix = [ upper(pVisual.pPlotError.sYScale) '_' ];
            else
                sPrefix = '';
            end
            sPrefix = [ sPrefix upper(pVisual.pPlotError.sLegendLocation) ];
            saveImage(hf, pVisual.pPlotError, sPrefix, mfilename); close(hf);
            if plotting('solution_norm') % plot norm of time-level solution
                Params = pVisual.pPlotError; 
                Params.sTitle    = {'Norm of solution'};
                Params.sLabelY   = {''};
                Params.sFileName = 'NORM_OF_SOLUTION';
                Params.nPosition = [0 0 0.8*768/1364 1];
                hf = plotErrorBehaviour(Params, ...
                    [ rError.CNormOfSolution(pCalc.nFreq); rError.L2NormOfSolution(pCalc.nFreq) ], []);
                saveImage(hf, Params, sPrefix, mfilename); close(hf);
            end
        end
    end
end
sFileName = [ pVisual.pPlotError.sPath '/' pVisual.pPlotError.sFilePrefix '_' pVisual.pPlotError.sFileName '_' pVisual.pPlotError.sFileSuffix ];
end

function [sLabelY, sFileSuffix] = prepareParams(pTask, pVisual, sFileSuffixE)
% function [sLabelY, sFileSuffix] = prepareParams(pTask, pCalc, pVisual, pCalc_E, sLabelYE, sFileSuffixE)
% sLabelYA = pVisual.pPlotError.sLabelY;
% if pCalc_E.n ~= pCalc.n || pCalc_E.m ~= pCalc.m
%     sLabelYE = [ sLabelYE ' (' ];
%     sLabelYA = [ sLabelYA ' (' ];
%     if pCalc_E.n ~= pCalc.n
%         sLabelYE = [ sLabelYE 'J=' num2str(pCalc_E.n-1) ];
%         sLabelYA = [ sLabelYA 'J=' num2str(pCalc.n  -1) ];
%     end
%     if pCalc_E.m ~= pCalc.m
%         if pCalc_E.n ~= pCalc.n
%             sLabelYE = [ sLabelYE ',' ];
%             sLabelYA = [ sLabelYA ',' ];
%         end
%         sLabelYE = [ sLabelYE 'M=' num2str(pCalc_E.m-1) ];
%         sLabelYA = [ sLabelYA 'M=' num2str(pCalc.m  -1) ];
%     end
%     sLabelYE = [ sLabelYE ')' ];
%     sLabelYA = [ sLabelYA ')' ];
% end
sLabelY = ''; % [ sLabelYE ' vs ' sLabelYA ]; %
if strcmp(sFileSuffixE, 'DTBC_VS_SDTBC') % strcmp(pVisual.pPlotError.sFileSuffix, sFileSuffixE)
    sFileSuffix = [ pVisual.pPlotError.sFileSuffix '_' sFileSuffixE ];
else
    sFileSuffix = [ pVisual.pPlotError.sFileSuffix '_VS_' sFileSuffixE ];
end
if isfield(pTask, 'V_Q')
    sFileSuffix = [ 'Q=' num2str(int32(max(pTask.V_Q))) '_' sFileSuffix ];
end
if isfield(pTask, 'V_lambda')
    sFileSuffix = [ 'LAMBDA=' num2str(pTask.V_lambda) '_' sFileSuffix ];
end
end