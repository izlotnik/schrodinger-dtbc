function Params = plot2DSolutionBehaviour(U, Params, pTask, k)

%% Initialize
if k == 0
    Params.sLegendTitle = {'$\left|\psi_{G}\right|$', 'Re $\psi_{G}$'};
    tx = ''; % [ 'Modulus and real part of \psi_{G} and normalized potential V' ]; %
    Params.sLabelY = '';
else
    tx = [ 'Modulus and real part of \Psi^m, m=' num2str(k) ];
end

%% Plot modulus and real part of U
% xlabel( Params.hAxes, tx );
% ylabel( Params.hAxes, Params.sLabelY );
plot(Params.hAxes, Params.nPlotX, abs(U) , Params.sLineType{1}...
                 , Params.nPlotX, real(U), Params.sLineType{2}...
                 , 'LineWidth', Params.nLineWidth );

%% Plot potential
if Params.bPlotPotential
    if strcmp(pTask.sExampleName, 'EX20r') || strcmp(pTask.sExampleName, 'EX22m') || strcmp(pTask.sExampleName, 'EX22r')
        xxx = Params.nPlotM;
    else
        xxx = Params.nPlotX;
    end
    P = V(pTask, xxx);
    if max(abs(P)) > 0
        P = P / max(abs(P));
        if strcmp(pTask.sExampleName, 'EX24')
            P = P * 0.5;
        end
    end
    set(Params.hAxes, 'NextPlot', 'add');
    plot(Params.hAxes, Params.nPlotX, P...
        , 'LineStyle', '-.', 'Color', [0 0.5 0], 'LineWidth', 2); % Params.sLineTypePotential - Params.nLineWidthPotential
    if strcmp(pTask.sExampleName, 'EX20') || strcmp(pTask.sExampleName, 'EX20r')
        Params.sLegendTitle(3) = {'$V_s$'};% {'$V/\max\left|V\right|$'}; 
    else
        Params.sLegendTitle(3) = {'$\chi_I$'};
    end
end
legend(Params.hAxes, Params.sLegendTitle, 'Location', Params.sLegendLocation, 'Interpreter', 'LaTeX');
set(Params.hAxes, 'NextPlot', Params.sNextPlot)
% chH = get(Params.hAxes, 'Children'); set(Params.hAxes, 'Children', [ chH(end); chH(1:end-1) ])
uistack( findobj(Params.hAxes, 'Color', 'b'), 'top');
% commented 2013/06/08 for video of EX22
% if max(max(abs(U))) > Params.nLimitY(2)
%     Params.nLimitY = 1.5 * Params.nLimitY;
% end
if isfield(Params, 'nLimitY')
    ylim(Params.nLimitY)
end

%% Pause
%%% EX20
% nBoundaryMultiplier = 1.2; % setCalcParams
% keyboard
% [hx,hy] = format_ticks(gca,{'$x_0=0$','','','','','','','','','$X_0$','$x_J=X$'},{'$-1$','','','','','$0$','','','','','$1$'});
% hold on; yL = get(gca, 'YLim'); xt = get(gca,'xtick'); line([xt(end) xt(end)], yL, '-k', 'LineWidth', Params.nLineWidthPotential);
if isfield(Params, 'bPlayPause') && Params.bPlayPause
    pause(1*0.01)
end

%% Movie control
Params.pMovie.nFrameCount = Params.pMovie.nFrameCount + 1;
Params.pMovie.pFrame(Params.pMovie.nFrameCount) = getframe(Params.hFigure);