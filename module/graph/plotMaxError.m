function sFileName = plotMaxError(e, pVisual, name, type, nM, nR, nJ, nN, iM, iR, iJ, iN)
p = plotInit;
sInd = type{2};
if strcmp(type{1}, 'rel'), sInd = [ sInd ',\,\rm ' type{1} ]; end
sFileSuffix_D = pVisual.sFileSuffix;
if iJ > 2 % Plot E^{N}(J) % iR == 2
    for m=1:(iM-1)
        hf = figure('units', pVisual.sUnits, 'outerposition', pVisual.nOuterPosition);
        set(gca, 'FontSize', pVisual.nFontSize, 'FontWeight', pVisual.sFontWeight);
        sFileSuffix = '';
        for N=1:(iN-1)
            [ sLegendTitle{N}, Nm, sFileSuffix, sFileSuffix_S ] = getParams(name, nN(N), N, sFileSuffix);
            yyy = e(m, 1, :, N); yyy = yyy(:);
            plot(1:size(e, 3), yyy ...
                , 'LineWidth',        pVisual.nLineWidth ...
                , 'LineStyle',        p.line{  mod(Nm-1,          length(p.line  ))+1} ...
                , 'Marker',           p.marker{mod(fix((Nm-1)/2), length(p.marker))+1} ...
                , 'MarkerFaceColor',  p.color( mod(Nm-1,          length(p.color ))+1, :) ...
                , 'Color',            p.color( mod(Nm-1,          length(p.color ))+1, :));
            set(gca, 'NextPlot', pVisual.sNextPlot);
        end
        legend(gca, sLegendTitle, 'Interpreter', 'LaTeX', 'Location', 'NorthEast'); % 'BestOutside'); % 'SouthWest'); %
        set(gca, 'YScale', 'log'); % set(gca, 'XScale', 'log');
        set(gca, 'FontSize', pVisual.nFontSize, 'FontWeight', pVisual.sFontWeight);
        % ylabel(['$\max\limits_{1\leq m\leq M} E_{' sInd '}^{m}$'], 'Interpreter', 'LaTeX');
        xlabel(pVisual.sLabelX, 'Interpreter', 'LaTeX');
        xlim([1 iJ-1]); % ylim([4*10^(-5) 1.5*10^(0)]); % ylim([5*10^(-6) 4*10^(-1)]);
        xt = get(gca, 'XTick');
        xti = fix(xt) == xt;
        xtn = NaN*ones(1, length(xt)); xtn(xti) = nJ(xt(xti));
        xtc = num2cell(xtn); xtc(isnan(xtn)) = {''};
        set(gca, 'XTickLabel', xtc);
        pVisual.sFileSuffix = [ type{1} '_' type{2} '_' sFileSuffix_S sFileSuffix(1:end-1) '_M=' num2str(nM(m)) '_' sFileSuffix_D ];
        sFileName = saveImage(hf, pVisual, [], mfilename); close(hf);
    end
else % Plot E^{R}(M) % iJ == 2, iN == 2
    hf = figure('units', pVisual.sUnits, 'outerposition', pVisual.nOuterPosition);
    set(gca, 'FontSize', pVisual.nFontSize, 'FontWeight', pVisual.sFontWeight);
    sFileSuffix = '';
    for R=1:(iR-1)
        [ sLegendTitle{R}, Nm, sFileSuffix, sFileSuffix_S ] = getParams('RICH', nR(R), R, sFileSuffix);
        plot(nM, e(:, R, 1, 1) ...
            , 'LineWidth',        pVisual.nLineWidth ...
            , 'LineStyle',        p.line{  mod(Nm-1, length(p.line  ))+1} ...
            , 'Marker',           p.marker{mod(Nm-1, length(p.marker))+1} ...
            , 'MarkerFaceColor',  p.color( mod(Nm-1, length(p.color ))+1, :) ...
            , 'Color',            p.color( mod(Nm-1, length(p.color ))+1, :));
        set(gca, 'NextPlot', pVisual.sNextPlot);
    end
    legend(gca, sLegendTitle, 'Interpreter', 'LaTeX', 'Location', 'NorthEast'); % 'BestOutside'); % 'SouthWest'); %
    set(gca, 'FontSize', pVisual.nFontSize, 'FontWeight', pVisual.sFontWeight);
    % ylabel(['$\max\limits_{1\leq m\leq M} E_{' sInd '}^{m}$'], 'Interpreter', 'LaTeX'); 
    set(gca, 'YScale', 'log'); set(gca, 'XScale', 'log');
    xlabel('M', 'Interpreter', 'LaTeX');
    xlim([min(nM) max(nM)]);
%     xlim([1 iM-1]);
%     xt = get(gca, 'XTick');
%     xti = fix(xt) == xt;
%     xtn = NaN*ones(1, length(xt)); xtn(xti) = nM(xt(xti));
%     xtc = num2cell(xtn); xtc(isnan(xtn)) = {''};
%     set(gca, 'XTickLabel', xtc);
    pVisual.sFileSuffix = [ type{1} '_' type{2} '_' sFileSuffix_S sFileSuffix(1:end-1) '_N=' num2str(nN(1)) '_J=' num2str(nJ(1)) '_' sFileSuffix_D ];
    sFileName = saveImage(hf, pVisual, [], mfilename); close(hf);
end
if iM > 2 && iR == 2 && iJ > 2 && iN == 2
    % close 'all'
    % set(0, 'defaulttextinterpreter', 'latex')
    Z = reshape(e, iM-1, iJ-1); Z(Z>=1) = 1.0;
    Contours = 10.^(-10:0);
    CC = {'10^{-10}','10^{-9}','10^{-8}','10^{-7}','10^{-6}','10^{-5}','10^{-4}','10^{-3}','10^{-2}','10^{-1}','10^{0}'};
    hf = figure('Units', pVisual.sUnits, 'OuterPosition', pVisual.nOuterPosition);
    set(gca, 'FontSize', pVisual.nFontSize, 'FontWeight', pVisual.sFontWeight);
    [X, Y] = meshgrid(nJ, nM);
    surf(X, Y, log10(Z)); box on; view(2); % contourf(X, Y, log10(Z), log10(Contours)); % , 'EdgeColor', 'None' % surf(X, Y, Z, Contours); box on; view(2); set(gca, 'ZScale', 'Log'); % grid on; % contourf(X, Y, log10(transpose(f)), log10(Contours), 'EdgeColor', 'None'); %
    set(gca, 'XTick', nJ, 'FontSize', pVisual.nFontSize, 'FontWeight', pVisual.sFontWeight); xlim([min(nJ) max(nJ)]); xlabel('J', 'FontSize', pVisual.nFontSize, 'FontWeight', pVisual.sFontWeight); 
    set(gca, 'YTick', nM, 'FontSize', pVisual.nFontSize, 'FontWeight', pVisual.sFontWeight); ylim([min(nM) max(nM)]); ylabel('M', 'FontSize', pVisual.nFontSize, 'FontWeight', pVisual.sFontWeight); 
    % surf(X, Y, log10(transpose(f))); % colorbar('Location', 'SouthOutside') % 'XTickLabel', nVal,
    colorbar('FontSize', pVisual.nFontSize, 'FontWeight', pVisual.sFontWeight, 'Location', 'SouthOutside', 'YTick', log10(Contours), 'YTickLabel', CC); % , 'XScale', 'Log' , 'XTick', log10(Contours), 'XTickLabel', CC
    % axis square % colormap(jet);
    caxis(log10([Contours(1) Contours(end)])); % set(c, 'YScale', 'Log'); % 
    % set(c, 'XTickLabel', 10.^get(c, 'XTick')) % set(c, 'XScale', 'Log')
    pVisual.sFileSuffix = [ type{1} '_' type{2} '_' 'CONTOUR' ];
    sFileName = saveImage(hf, pVisual, [], mfilename); close(hf);
end
end

function [ sLegendTitle, Nm, sFileSuffix, sFileSuffix_S ] = getParams(name, N, iN, sFileSuffix)
switch name
    case 'RICH'
        sLegendTitle  = [ '$r='       num2str(N) '$' ];
        sFileSuffix   = [ sFileSuffix num2str(N) ',' ];
        sFileSuffix_S = 'R=';
        Nm = N + 10;
    case 'FEM'
        sLegendTitle  = [ '$n='       num2str(N) '$' ];
        sFileSuffix   = [ sFileSuffix num2str(N) ',' ];
        sFileSuffix_S = 'N=';
        Nm = N;
    case 'FFDS'
        sLegendTitle  = [ '$\theta='         sym2str(sym(N)) '$'          ];
        sFileSuffix   = [ sFileSuffix strrep(sym2str(sym(N)),'/','-') ',' ];
        sFileSuffix_S = 'PARAM=';
        switch iN
            case 1/4
                Nm = 17;
            case 1/6
                Nm = 1; % ~ n=1
            case 1/12
                Nm = 3; % ~ n=3
            case 0
                Nm = 12;
                % theta < 0
            case -1/12
                Nm = 13;
            case -1/6
                Nm = 14;
            case -1/4
                Nm = 15;
            case -1/2
                Nm = 16;
                % theta > 1/4
            case 1/2
                Nm = 18;
            case 1
                Nm = 19;
            case 2
                Nm = 20;
            otherwise
                Nm = N;
        end
end
end