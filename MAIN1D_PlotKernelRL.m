clearContext; format short; echo on; tic
%%
sExampleName   = 'EX01';                % 'EX02m'; %
sDataDirectory = 'C:\MATLAB\R2011b';    % 'C:\Matlab\R2009b'; %  'E:\MATLAB\R2011b'; %
sBaseRunID     = '';                    % ''; % 'LAST'; % '20100205_170006_125'; %
[pTask pInitWave] = setTaskParams('Schrodinger', '1D', sExampleName, sDataDirectory, sBaseRunID);
%% EXO1
nM = 3000; nM_max = 1000; % nM;
nJ = 30;
NN = 1:9; % {1:6,[6 7],[6 8],[6 9]};
sT = 'FEM';
%% TEST
% nM = 350; nM_max = nM; % 1000; % 350; %
% nJ = 55;
% NN = 7; % [5 7]; % [1/6 1/12 0 1/4]; % 1:10
%%
if strcmp(sT, 'FFDS')
    sMethod = {'TRISYS'};
else
    sMethod = {'QR'};
end
K = NaN * ones(length(1:nM), length(NN));
%%
N1 = 1;
for N=NN
    [pCalc pFile] = setCalcParams(pTask, {sT}, [N nJ nM NaN 0], sMethod);
    [pVisual pCalc] = setVisualParams(pTask, pCalc, [0 0 0 0 0], [0 0 0 0], [1 1]);
    N2 = fix((N+1)/2); % R1 = NaN * ones(nM, N2); L1 = R1;
    pTask.V_inf = pTask.V_Rinf;
    [ K1, c_0, R1, L1 ] = calcConvolutionKernel(pTask, pCalc);
    %%
    p = plotInit;
    pVisual.pPlotKernel.sLabelY  = ['$\left|R^{(' num2str(NN(N1)) '),m}_{\ell}\right|$'];
    hf = figure('units', pVisual.pPlotKernel.sUnits, 'outerposition', pVisual.pPlotKernel.nOuterPosition);
    sFileSuffix = 'R';
    for N3=1:N2
        sLegendTitle{N3} = ['$\ell=' num2str(N3) '$'];
        plot(1:min(nM,nM_max), abs(R1(1:min(nM,nM_max),N3))...
            ,'LineWidth',       pVisual.pPlotKernel.nLineWidth ...
            ,'LineStyle',       p.line{  mod(N3-1,length(p.line  )-1)+1} ...
            ,'Color',           p.color( mod(N3-1,length(p.color ))+1,:));
        set(gca, 'NextPlot', pVisual.pPlotKernel.sNextPlot);
    end
    legend(gca, sLegendTitle, 'Location', pVisual.pPlotKernel.sLegendLocation, 'Interpreter', 'LaTeX');
    set(gca, 'YScale', 'log'); % set(gca, 'XScale', 'log');
    set(gca, 'FontSize', pVisual.pPlotKernel.nFontSize, 'FontWeight', pVisual.pPlotKernel.sFontWeight);
    xlabel(pVisual.pPlotKernel.sLabelX, 'Interpreter', 'LaTeX');
    ylabel(pVisual.pPlotKernel.sLabelY, 'Interpreter', 'LaTeX');
    switch pCalc.name
        case 'FEM'
            pVisual.pPlotKernel.sFileSuffix = [sFileSuffix '_N='     num2str(NN(N1))                     ];
        case 'FFDS'
            pVisual.pPlotKernel.sFileSuffix = [sFileSuffix '_PARAM=' strrep(sym2str(sym(NN(N1))),'/','-')];
    end
    saveImage(hf, pVisual.pPlotKernel, ['M=' num2str(nM_max) '_' 'J=' num2str(nJ) '_' 'MMAX=' num2str(min(nM,nM_max))], mfilename); % close(hf);
    clear sLegendTitle;
    %%
    if mod(NN(N1),2) ~= 0 % N - нечетное
        N2 = N2 - 1;
    end
    if N2 > 0
        pVisual.pPlotKernel.sLabelY  = ['$\left|L^{(' num2str(NN(N1)) '),m}_{\ell}\right|$'];
        hf = figure('units', pVisual.pPlotKernel.sUnits, 'outerposition', pVisual.pPlotKernel.nOuterPosition);
        sFileSuffix = 'L';
        for N3=1:N2
            sLegendTitle{N3} = ['$\ell=' num2str(N3) '$'];
            plot(1:min(nM,nM_max),abs(L1(1:min(nM,nM_max),N3))...
                ,'LineWidth',       pVisual.pPlotKernel.nLineWidth ...
                ,'LineStyle',       p.line{  mod(N3-1,length(p.line  )-1)+1} ...
                ,'Color',           p.color( mod(N3-1,length(p.color ))+1,:));
            set(gca, 'NextPlot', pVisual.pPlotKernel.sNextPlot);
        end
        legend( gca, sLegendTitle, 'Location', pVisual.pPlotKernel.sLegendLocation, 'Interpreter', 'LaTeX');
        set(gca, 'YScale', 'log'); % set(gca, 'XScale', 'log');
        set(gca, 'FontSize', pVisual.pPlotKernel.nFontSize, 'FontWeight', pVisual.pPlotKernel.sFontWeight);
        xlabel(pVisual.pPlotKernel.sLabelX, 'Interpreter', 'LaTeX');
        ylabel(pVisual.pPlotKernel.sLabelY, 'Interpreter', 'LaTeX');
        switch pCalc.name
            case 'FEM'
                pVisual.pPlotKernel.sFileSuffix = [sFileSuffix '_N='     num2str(NN(N1))                     ];
            case 'FFDS'
                pVisual.pPlotKernel.sFileSuffix = [sFileSuffix '_PARAM=' strrep(sym2str(sym(NN(N1))),'/','-')];
        end
        saveImage(hf, pVisual.pPlotKernel, ['M=' num2str(nM_max) '_' 'J=' num2str(nJ) '_' 'MMAX=' num2str(min(nM,nM_max))], mfilename); % close(hf);
        clear sLegendTitle;
    end
    N1 = N1 + 1;
end
toc