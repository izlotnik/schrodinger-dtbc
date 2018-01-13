clearContext; format short; echo on; tic
%%
sExampleName   = 'EX01';                % 'EX02m'; %
sDataDirectory = 'C:\MATLAB\R2011b';    % 'C:\Matlab\R2009b';    % 'D:\MATLAB\R2011b'; %
sBaseRunID     = '';                    % ''; % 'LAST'; % '20100205_170006_125'; %
[pTask pInitWave] = setTaskParams('Schrodinger', '1D', sExampleName, sDataDirectory, sBaseRunID);
%% FEM
nM = 3000; nM_max = 1000; % nM; %
nJ = 30;
NN = {1:6,[6 7],[6 8],[6 9]};
sT = 'FEM';
%% EX01 (FFDS) A
% nM = 300; nM_max = nM; %
% nJ = 240;
% NN = {[1/4 1/6 1/12 0]};
% sT = 'FFDS';
%% EX01 (FFDS) B
% nM = 3000; nM_max = 50; % nM;
% nJ = 1000; % 30; %
% NN = {[1/4 1/6 1/12 0]};
% sT = 'FFDS';
%%
if strcmp(sT, 'FFDS')
    sMethod = {'TRISYS'};
else
    sMethod = {'QR'};
end
K = NaN * ones(nM, length([NN{:}]));
N1 = 1;
for N=[NN{:}]
    [pCalc pFile] = setCalcParams(pTask, {sT}, [N nJ nM NaN 0], sMethod);
    [pVisual pCalc] = setVisualParams(pTask, pCalc, [0 0 0 0 0], [0 0 0 0], [1 1]);
    pTask.V_inf = pTask.V_Rinf;
    K1 = calcConvolutionKernel(pTask, pCalc);
    K(:, N1)= K1(1:nM);
    N1 = N1 + 1;
end
toc
nM_Default = nM; nM = min(nM, nM_max);
%%
p  = plotInit;
N1 = 1;
for q=1:length(NN)
    nNc = NN{q};
    hf = figure('units', pVisual.pPlotKernel.sUnits, 'outerposition', pVisual.pPlotKernel.nOuterPosition);
    sFileSuffix = '';
    for N=1:length(nNc)%1:(N1-5)%
        switch pCalc.name
            case 'FEM'
                sLegendTitle{N} = [ '$n='       num2str(nNc(N)) '$' ];
                sFileSuffix     = [ sFileSuffix num2str(nNc(N)) ',' ];
                Nm = nNc(N);
            case 'FFDS'
                sLegendTitle{N} = [ '$\theta='         sym2str(sym(nNc(N)))            '$' ];
                sFileSuffix     = [ sFileSuffix strrep(sym2str(sym(nNc(N))), '/', '-') ',' ];
                switch nNc(N)
                    case 1/4,   Nm = 17;
                    case 1/6,   Nm =  1; % ~ n=1
                    case 1/12,  Nm =  3; % ~ n=3
                    case 0,     Nm = 12; % theta < 0
                    case -1/12, Nm = 13;
                    case -1/6,  Nm = 14;
                    case -1/4,  Nm = 15;
                    case -1/2,  Nm = 16; 
                    case 1/2,   Nm = 18; % theta > 1/4
                    case 1,     Nm = 19;
                    case 2,     Nm = 20;
                    otherwise,  Nm = nNc(N);
                end
        end
        if Nm==17 % theta = 1/4
            K(2:2:nM, N1) = K(1:2:nM, N1);
        end
        plot(1:nM, abs(K(1:nM, N1))...
            ,'LineWidth',       pVisual.pPlotKernel.nLineWidth ...
            ,'LineStyle',       p.line{  mod(Nm-1, length(p.line  ))+1} ...
            ,'Color',           p.color( mod(Nm-1, length(p.color ))+1, :));
        set(gca, 'NextPlot', pVisual.pPlotKernel.sNextPlot);
        N1 = N1 + 1;
    end
    legend(gca, sLegendTitle, 'Location', pVisual.pPlotKernel.sLegendLocation, 'Interpreter', 'LaTeX');
    set(gca, 'YScale', 'log');
    set(gca, 'FontSize', pVisual.pPlotKernel.nFontSize, 'FontWeight', pVisual.pPlotKernel.sFontWeight);
    xlabel('$m$', 'Interpreter', 'LaTeX');
    ylabel('$\left| K^{(n),m}_{\rm ref}\right|$', 'Interpreter', 'LaTeX');
    switch pCalc.name
        case 'FEM'
            pVisual.pPlotKernel.sFileSuffix = [ 'N='     sFileSuffix(1:end-1) ];
        case 'FFDS'
            pVisual.pPlotKernel.sFileSuffix = [ 'PARAM=' sFileSuffix(1:end-1) ];
    end
    saveImage(hf, pVisual.pPlotKernel, [ 'M=' num2str(nM_Default) '_' 'J=' num2str(nJ) '_' 'MMAX=' num2str(nM_max) ], mfilename); % close(hf);
    clear sLegendTitle;
end