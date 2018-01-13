clearContext; format short; echo off;
%%
[pTask, pInitWave] = setTaskParams('Schrodinger', '1D', 'EX01r', 'D:\Matlab\R2011b', '');
% [pTask, pInitWave] = setTaskParams('Schrodinger', '1D', 'EX01r', 'E:\Matlab\R2011b', '');
%%
% conv      - 17'205 sec: NN = 1:10;      mm = round(10.^(1:0.5:5)); q_max = 5;
% conv      -  8'513 sec: NN = [1:5 7 9]; mm = round(10.^(1:0.5:5)); q_max = min(max(11-N, 5), 10);
% conv_fft2 -    173 (PC-NIX:116, TESLA:57) sec: NN = [1:5 7 9]; mm = round(10.^(1:0.5:5)); q_max = min(max(11-N, 5), 10);
n  = 30; sConvMethod = 'conv2'; % 'conv\_fft2'; % 
NN = [1 3 5 7 9];           % 1;    % 1:10; % [1/6 1/12 0 1/4]; %
mm = round(10.^(2:5));  % 10^7; % round(10.^(1:0.5:5)); % mm = sort([mm mm(2:end)/2]);
t  = zeros(length(mm), length(NN));
m1 = 1;
for m=mm
    N1 = 1;
    for N=NN
        [ pCalc,   pFile ] = setCalcParams(pTask, {'FEM'}, [N n m NaN 0], {'QR'});
        [ pVisual, pCalc ] = setVisualParams(pTask, pCalc, [0 0 0 0 0], [0 0 0 0], [1 1]);
        pLTask = pTask; pLTask.V_inf = pTask.V_Linf;
        pLCalc = pCalc; pLCalc.name = pCalc.pLBoundary.sName;
        switch pLCalc.name
            case 'FFDS'
                pLCalc.theta = pCalc.pLBoundary.theta;
            case 'FEM'
                pLCalc.N = pCalc.pLBoundary.N;
        end
        q_max = 1; % min(max(11-N, 5), 10); % 15; % min(max(16-N, 7), 11); % 
        t_res = NaN * ones(1, q_max);
        for q=1:q_max
            tic
            [ K, c ] = calcConvolutionKernel(pLTask, pLCalc);
            t_res(q) = toc;
        end
        % t_res = sort(t_res); t_res = t_res(2:end-1);
        t(m1, N1) = mean(t_res);
        N1 = N1 + 1;
        clear pLTask pLCalc;
    end
    m1 = m1 + 1;
end
sum(sum(t)) * 5

%%
tic
clear sLegendTitle;
p  = plotInit;
hf = figure('units', pVisual.pPlotKernel.sUnits, 'outerposition', pVisual.pPlotKernel.nOuterPosition);
set(gca, 'FontSize', pVisual.pPlotKernel.nFontSize, 'FontWeight', pVisual.pPlotKernel.sFontWeight);
sFileSuffix = '';
pVisual.pPlotKernel.sLabelX         = 'm';
pVisual.pPlotKernel.sLabelY         = 'CPU time, sec.';
pVisual.pPlotKernel.sLegendLocation = 'NorthWest';
N1 = 1;
for N=NN % [1:5 7] % 
    switch pCalc.name
        case 'FEM'
            sLegendTitle{N1} = [ '$n='        num2str(NN(N1)) '$' ];
            sFileSuffix      = [ sFileSuffix  num2str(NN(N1)) ',' ];
        case 'FFDS'
            sLegendTitle{N1} = [ '$\theta='         sym2str(sym(NN(N1)))            '$' ];
            sFileSuffix      = [ sFileSuffix strrep(sym2str(sym(NN(N1))), '/', '-') ',' ];
    end
    loglog( mm, t(:, N1)... % mm(3:end), t(3:end, N)... % 
        ,'LineWidth',       pVisual.pPlotKernel.nLineWidth ...
        ,'LineStyle',       p.line{  mod(N-1, length(p.line  ))+1} ...
        ,'Color',           p.color( mod(N-1, length(p.color ))+1, :));
    set(gca, 'nextplot', pVisual.pPlotKernel.sNextPlot);
    N1 = N1 + 1;
end
set(gca, 'FontSize', pVisual.pPlotKernel.nFontSize, 'FontWeight', pVisual.pPlotKernel.sFontWeight); % set(gca, 'XTickLabel', sprintf('%2.3f|', nFreq(int32(get(gca, 'XTick')))));
legend(gca, sLegendTitle, 'Location', pVisual.pPlotKernel.sLegendLocation, 'Interpreter', 'LaTeX');
% xlabel('M'); ylabel(pVisual.pPlotKernel.sLabelY);
xlim([ min(mm) max(mm) ]); % ylim([ 10^(-4) 10^(3) ]);
% title('The time of computation for DTBC kernel $K^{(n),m}_{\rm ref}$', 'Interpreter', 'LaTeX');
switch pCalc.name
    case 'FEM'
        sFileSuffix = [ 'N='     sFileSuffix(1:end-1) ];
    case 'FFDS'
        sFileSuffix = [ 'PARAM=' sFileSuffix(1:end-1) ];
end
pVisual.pPlotKernel.sFileSuffix = [ strrep(upper(sConvMethod), '\', '') '_' sFileSuffix ];
saveImage(hf, pVisual.pPlotKernel, [ 'M=' num2str(max(mm)) ], mfilename); % close(hf);
csvwrite([ pVisual.pPlotKernel.sPath '/' pVisual.pPlotKernel.sFilePrefix '_' pVisual.pPlotKernel.sFileName '_' pVisual.pPlotKernel.sFileSuffix '.csv' ], t);
toc