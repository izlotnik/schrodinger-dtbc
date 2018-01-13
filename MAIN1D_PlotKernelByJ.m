clearContext; format short; echo on; tic
%%
sExampleName   = 'EX01';                % 'EX02m'; %
sDataDirectory = 'C:\Matlab\R2009b';    % 'E:\MATLAB\R2011b'; % 'D:\MATLAB\R2011b'; % 
sBaseRunID     = '';                    % ''; % 'LAST'; % '20100205_170006_125'; %
[pTask pInitWave] = setTaskParams( 'Schrodinger', '1D', sExampleName, sDataDirectory, sBaseRunID );
%%
m = 5000; mMax = m; % 1000;
nn = 10:30:100;
N = 1/12; % 1/4; % 
K = NaN * ones(m,length(nn));
n1 = 1;
for n=nn
    %%
    [pCalc pFile] = setCalcParams(pTask, {'FFDS'}, [N n m NaN 0], {'QR'});
    [pVisual pCalc] = setVisualParams(pTask, pCalc, [0 0 0 0 0], [0 0 0 0], [1 1]);
    pTask.V_inf = pTask.V_Rinf;
    K1 = calcConvolutionKernel(pTask, pCalc);
    K(:,n1)= K1(1:m);
    n1 = n1 + 1;
end
[pCalc pFile] = setCalcParams(pTask, {'FEM','SD'}, [N n m NaN 0], {'QR'});
[pVisual pCalc] = setVisualParams(pTask, pCalc, [0 0 0 0 0], [0 0 0 0], [1 1]);
pRTask = pTask; pRCalc = pCalc;
pRTask.V_inf = pTask.V_Rinf;
pRCalc.name = pCalc.pRBoundary.sName;
switch pRCalc.name
    case 'FFDS'
        pRCalc.theta = pCalc.pRBoundary.theta;
    case 'FEM'
        pRCalc.N = pCalc.pRBoundary.N;
end
K1 = calcConvolutionKernel(pRTask, pRCalc);
K(:,n1)= K1(1:m);
toc
%%
p = plotInit;
hf = figure('units', pVisual.pPlotKernel.sUnits, 'outerposition', pVisual.pPlotKernel.nOuterPosition);
switch pCalc.name
    case 'FEM'
        sLabelY = ['$n=' sym2str(sym(N)) '$'];
        sFileSuffix = num2str(N);
    case 'FFDS'
        sLabelY = ['$\theta=' sym2str(sym(N)) '$'];
        sFileSuffix = strrep(sym2str(sym(N)),'/','-');
end
for n=1:(n1-1)
    sLegendTitle{n} = ['$J=' num2str(nn(n)) '$'];
    plot(1:min(m,mMax), abs(K(1:min(m,mMax),n))...
            ,'LineWidth',       pVisual.pPlotKernel.nLineWidth ...
            ,'LineStyle',       p.line{  mod(n-1,length(p.line  ))+1} ...
            ,'Color',           p.color( mod(n-1,length(p.color ))+1,:));
    set(gca, 'nextplot', pVisual.pPlotKernel.sNextPlot);
end
sLegendTitle{n1} = '$SD$';
plot(1:min(m,mMax),abs(K(1:min(m,mMax),n1))...
    ,'LineWidth',       10 ...
    ,'LineStyle',       p.line{  mod(n1-1,length(p.line  ))+1} ...
    ,'Color',           p.color( mod(n1-1,length(p.color ))+1,:));
legend(gca, sLegendTitle, 'Location', pVisual.pPlotKernel.sLegendLocation, 'Interpreter', 'LaTeX');
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');
set(gca, 'FontSize', pVisual.pPlotKernel.nFontSize, 'FontWeight', pVisual.pPlotKernel.sFontWeight);
xlabel(pVisual.pPlotKernel.sLabelX);
ylabel(sLabelY, 'Interpreter', 'LaTeX');
saveImage( hf, pVisual.pPlotKernel, ['M=' num2str(min(m,mMax))], mfilename ); % close(hf);
clear sLegendTitle;