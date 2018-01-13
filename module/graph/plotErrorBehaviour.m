function [hf, ax] = plotErrorBehaviour(Params, vAbsErrors, vRelErrors)
hf = figure('Name', Params.sName...
            , 'NumberTitle'     , Params.sNumberTitle...
            , 'Units'           , Params.sUnits...
            , 'OuterPosition'   , Params.nPosition);
N = ~isempty(vAbsErrors)+~isempty(vRelErrors);
ax = NaN * ones(1, N);
ylim_c = [];

m = 23;
% subplot(m, 2, 1:2:(2*m-3))
% l=plot(0:.01:4,sin(0:.01:4));
% subplot(m, 2, 2:2:(2*m-2))
% r=plot(0:.01:4,cos(0:.01:4),'r');
% A = subplot(m, 2, (2*m-1):(2*m));
% legend(A,[l r], 'sin', 'cos', 'location', 'best')
% set(A, 'visible', 'off')

for i=1:N
    if i==1
        ppp = 1:((m-1)/2-3);
    else
        ppp = ((m+1)/2+1):(m-3);
    end
    ax(i) = subplot(1, m, ppp... % 1, N, i... % m, 2, i:2:(2*m-4+i)... % 
        , 'FontSize'       , Params.nFontSize   ...
        , 'FontWeight'     , Params.sFontWeight ...
        , 'NextPlot'       , 'Add');
    switch i
        case 1, data = vAbsErrors;
        case 2, data = vRelErrors;
    end
    for j=1:size(data, 1)
        if isfield(Params, 'sLineStyle') && isfield(Params, 'sColor')
            p(i, j) = plot(data(j, :), 'Parent', ax(i) ...
                , 'LineWidth',        Params.nLineWidth ...
                , 'LineStyle',        Params.sLineStyle{j} ...
                , 'Marker',           'n' ...
                , 'MarkerFaceColor',  Params.sColor{j} ...
                , 'Color',            Params.sColor{j});
        else
            p(i, j) = plot(data(j, :), Params.sLineType{j}, 'Parent', ax(i), 'LineWidth', Params.nLineWidth);
        end
    end
    xlim(ax(i), Params.nLimitX);
    lll = length(data(j, :))-1;
    xt = 0:(lll/4):lll; % 0:(lll/3):lll; %
    set(ax(i), 'XTick', xt); set(ax(i), 'XTickLabel', sprintf('%.3g|', Params.nFreq(end)*xt/(Params.nLimitX(end)-1)));
    % set(ax(i), 'XTickLabel', sprintf('%2.3f|', Params.nFreq(1+int32(get(gca, 'XTick')))));
    % set(ax(i), 'XTickLabel', Params.nFreq(1+int32(get(gca, 'XTick'))), 'FontSize', Params.nFontSize-4);
    if isfield(Params, 'sLabelY') && i == 1
        ylabel(ax(i), Params.sLabelY); 
    end
%     if isfield(Params, 'sYScale')
%         set(ax(i), 'YScale', Params.sYScale);
%     end
    set(ax(i), 'YScale', 'log');
    if isfield(Params, 'nLimY')
        set(ax(i), 'ylim', Params.nLimY)
        yt = Params.nLimY(1)*10.^(0:log10(Params.nLimY(2)/Params.nLimY(1)));
        for j = 1:length(yt)
            ytl{j} = [ '10^{' num2str(log10(yt(j))) '}' ];
        end
        set(ax(i), 'YTick', yt); % set(ax(i), 'YTickLabel', ytl);
        % tick2text(ax(i), 'axis', 'y', 'yformat', @(x) sprintf('10^{%d}', log10(x)));
        % set(getappdata(ax(i), 'YTickText'), 'FontSize', Params.nFontSize, 'FontWeight', Params.sFontWeight);
    end
    if isfield(Params, 'sTitle') && ~isempty(Params.sTitle{1})
        title(ax(i), Params.sTitle{i});
    end
    ylim_c = [ ylim_c get(ax(i), 'ylim') ];
end
if Params.bAlignYLim
    set(ax, 'ylim', [ min(ylim_c) max(ylim_c) ])
end

%% Add legend
A = subplot(1, m, (m-2):m); % A = subplot(m, 2, (2*m-1):(2*m));
legend(A, p(1, :), Params.sLegendTitle, 'Location', Params.sLegendLocation, 'Interpreter', 'LaTeX', 'FontSize', Params.nFontSize+8, 'FontWeight', Params.sFontWeight);
set(A, 'Visible', 'Off')
% legend(ax(2), Params.sLegendTitle, 'Location', Params.sLegendLocation, 'Interpreter', 'LaTeX');
% for i=1:N
%     legend(ax(i), Params.sLegendTitle, 'Location', Params.sLegendLocation, 'Interpreter', 'LaTeX'); % , 'Color', 'none', 'Box', 'off'
% end
end