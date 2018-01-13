function plotKernel(R,pCalc,pVisual)
%%
mMax = 500;
mm = 1:(pCalc.m-1); % mm = 1:((pCalc.m-1)/25);
hf = figure('units', pVisual.sUnits, 'outerposition', pVisual.nOuterPosition);
set(gca, 'FontSize', pVisual.nFontSize, 'FontWeight', pVisual.sFontWeight);
% nMultiplier = 1.5;
% hhh = plot(log10(mm),log10(abs(R(mm))),'-', log10(mm), log10(nMultiplier*mm.^(-3/2)),'--');
hhh = plot(mm,log10(abs(R(mm)))); % plot(log10(mm),log10(abs(R(mm))));
% hhh = plot(mm, abs(R(mm)),'-', mm, nMultiplier*mm.^(-3/2),'--');
set(hhh, 'LineWidth', pVisual.nLineWidth);
legend(gca, pVisual.sLegendTitle, 'Interpreter', 'LaTeX');
% legend(gca, '|R^m|', 'Re R^m', 'Im R^m', '1.5 m^{-3/2}', 'Interpreter', 'LaTeX');
xLim = [min(mm) min(max(mm),mMax)]; % round([min(log10(mm)) max(log10(mm))]*10)/10; % [0.25 2.55]; %
xlim(xLim);
ylabel(pVisual.sLabelY);
xlabel('m'); % xlabel([num2str(min(mm)-1) '\leq m \leq ' num2str(min(max(mm),mMax))]);
saveImage( hf, pVisual, ['M=' num2str(min(max(mm),mMax))], mfilename )
% close(hf)