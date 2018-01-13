clearContext; format short; echo on;
%%
sExampleName   = 'EX22a';               % 'EX01'; % 'EX23'; % 
sDataDirectory = 'C:\MATLAB\R2011b';    % 'C:\Matlab\R2009b';    % 'D:\MATLAB\R2011b'; % 
sBaseRunID     = 'LAST';                % ''; % 'LAST'; % '20100205_170006_125'; %
[ pTask pInitWave ] = setTaskParams('Schrodinger', '2D', sExampleName, sDataDirectory, sBaseRunID); pTaskE = pTask;
%% Calc "exact" solution
[pCalc   pFile] = setCalcParams(pTask, {'FFDS'}, [0 4*1200/1 4*2*(3/5)*1000/1 NaN 0 0 4*2*64/1], {'TRISYS', 'SSP', 'CONST'}); % 
% [pCalc   pFile] = setCalcParams(pTask, {'FFDS'}, [0 1*1200/4/10 1*(3/5)*1000/2 NaN 0 0 1*128/4], {'TRISYS', 'SSP'}); 
% EX02m : 4*1200/1 4*(3/5)* 800/1 NaN 0 0 4*128/1
% EX22a : 8*1200/1 8*(3/5)*1000/1 NaN 0 0 4*128/1
[pVisual pCalc] = setVisualParams(pTask, pCalc, [0 0 0 0 0], [1 10 0 0], [0 1]); % [1 0 10 0 0]
pCalc.bCalcError = false; % pVisual.bPlotTimeBehaviour = false;
echo off; tic
MAIN2D_ONE_TIME_LEVEL
echo on ; toc
pCalcE = pCalc; pFileE = pFile; sLabelYE_D = pVisual.pPlotError.sLabelY; sFileSuffixE = pVisual.pPlotError.sFileSuffix;
P = [ ...
%        4*1200/1 4*(3/5)*1000/1 4*64/1   
%        2*1200/1 2*(3/5)*1000/1 2*64/1
        1*1200/1 2*(3/5)*1000/1 2*64/1
%        1*1200/1 1*(3/5)*1000/1 1*64/1
%        1*1200/4/10 1*(3/5)*1000/2 1*128/4
%        8*1200/1 8*(3/5)*1000/1 8*64/1
%        1*1200/4 1*(3/5)*1000/4 1*128/4
%        1*1200/1 1*(3/5)*1000/1 1*128/2
%        1*1200/1 2*(3/5)*1000/1 1*128/1       
%        1*1200/1 4*(3/5)*1000/1 1*128/1       
    ];
% s = size(P, 1)+1;
% for j=2.^(-4:1)%(-1:1)
%     for m=2.^0%(-1:1)
%         for k=2.^(-4:1)%(-1:1)
%             P(s, :) = [ j*1200/1 m*(3/5)*1000/1 k*64/1 ];
%             s = s + 1;
%         end
%     end
% end
%%
% s = size(P, 1)+1;
% for j=2.^(-4:-3)%(-2:2)
%     for m=2.^(2)
%         for k=2.^(2)
%             P(s, :) = [ j*1200/1 m*(3/5)*1000/1 k*64/1 ];
%             s = s + 1;
%         end
%     end
% end
% s = size(P, 1)+1;
% for j=2.^(2)
%     for m=2.^(-3:-3)%(-2:2)
%         for k=2.^(2)
%             P(s, :) = [ j*1200/1 m*(3/5)*1000/1 k*64/1 ];
%             s = s + 1;
%         end
%     end
% end
% s = size(P, 1)+1;
% for j=2.^(2)
%     for m=2.^(2)
%         for k=2.^(-4:-3)%(-2:2)
%             P(s, :) = [ j*1200/1 m*(3/5)*1000/1 k*64/1 ];
%             s = s + 1;
%         end
%     end
% end
% P = unique(P, 'rows');
for pp=1:size(P, 1)
    pTask = pTaskE; pCalc = []; pFile = []; pVisual = []; pp
    % sExampleName  = 'EX22a'; [ pTask pInitWave ] = setTaskParams('Schrodinger', '2D', sExampleName, sDataDirectory, sBaseRunID); pTaskE = pTask;
    %% Calc solution
    [pCalc   pFile] = setCalcParams(pTask, {'FFDS'}, [0 P(pp, 1) P(pp, 2) NaN 0 0 P(pp, 3)], {'TRISYS', 'SSP', 'CONST'}); % 
    [pVisual pCalc] = setVisualParams(pTask, pCalc, [1 1 10 0 0], [1 10 0 0], [0 1]); % [1 0 10 0 0], [1 10 0 0]
    pCalc.bCalcError = false; % pVisual.bPlotTimeBehaviour = false;
    echo off; tic
    MAIN2D_ONE_TIME_LEVEL
    P(pp, 8) = toc;
    if pVisual.bPlotTimeBehaviour
        close(pVisual.pPlot3D.hFigure);
    end
    pCalcA = pCalc; pFileA = pFile; sLabelYA = pVisual.pPlotError.sLabelY; sFileSuffixA = pVisual.pPlotError.sFileSuffix;
    %% Calc absolute and relative errors
    tic
    sLabelYE = sLabelYE_D;
    if pCalcE.n ~= pCalcA.n || pCalcE.K ~= pCalcA.K || pCalcE.m ~= pCalcA.m
        sLabelYE = [ sLabelYE ' (' ];
        sLabelYA = [ sLabelYA ' (' ];
        if pCalcE.n ~= pCalcA.n
            sLabelYE = [ sLabelYE 'J=' num2str(pCalcE.n-1) ];
            sLabelYA = [ sLabelYA 'J=' num2str(pCalcA.n-1) ];
        end
        if pCalcE.K ~= pCalcA.K
            if pCalcE.n ~= pCalcA.n
                sLabelYE = [ sLabelYE ',' ];
                sLabelYA = [ sLabelYA ',' ];
            end
            sLabelYE = [ sLabelYE 'K=' num2str(pCalcE.K-1) ];
            sLabelYA = [ sLabelYA 'K=' num2str(pCalcA.K-1) ];
        end
        if pCalcE.m ~= pCalcA.m
            if pCalcE.K ~= pCalcA.K || pCalcE.n ~= pCalcA.n
                sLabelYE = [ sLabelYE ',' ];
                sLabelYA = [ sLabelYA ',' ];
            end
            sLabelYE = [ sLabelYE 'M=' num2str(pCalcE.m-1) ];
            sLabelYA = [ sLabelYA 'M=' num2str(pCalcA.m-1) ];
        end
        sLabelYE = [ sLabelYE ')' ];
        sLabelYA = [ sLabelYA ')' ];
    end
    s1 = (pCalcE.K-1)/(pCalcA.K-1);
    s2 = (pCalcE.n-1)/(pCalcA.n-1);
    s3 = (pCalcE.m-1)/(pCalcA.m-1);
    h = pCalc.h_mod; delta = pCalc.delta;
    er1 = zeros(1, pCalcA.m); er1_e = er1; er2 = er1; er2_e = er1;
    hw = waitbar(0, 'Calc error: please wait...');
    for k=pCalc.nFreq
        if k > 1
            sFullName = [ pCalcE.sWorkDirectory '\' pFileE.pSolution.sName '_M=' num2str((k-1)*s3) '.' pFileE.pSolution.sExtension ];
        else
            sFullName = pFileE.pSolution.sFullName;
        end
        U_E = loadVariable(sFullName, pFileE.pSolution.bLoad, mfilename); U_E = [ zeros(1, pCalcE.n); U_E; zeros(1, pCalcE.n) ];
        if k > 1
            sFullName = [ pCalcA.sWorkDirectory '\' pFileA.pSolution.sName '_M=' num2str(k-1) '.' pFileA.pSolution.sExtension ];
        else
            sFullName = pFileA.pSolution.sFullName;
        end
        U_A = loadVariable(sFullName, pFileA.pSolution.bLoad, mfilename); U_A = [ zeros(1, pCalcA.n); U_A; zeros(1, pCalcA.n) ];
        %%
%         if k==1
%             pVisual.pPlot3D = initGraphics('3D', pVisual.pPlot3D);
%         end
%         Ua = U_E(1:s1:end, 1:s2:end) - U_A(:, :); 
%         % x = pCalc.h_mod*((0:(pCalc.n_mod-1))); y = pCalc.delta*(1:(pCalc.K-2));
%         % ind_x = 1.5-0.5 <= x & x <= 1.6+0.5;
%         % Ua(:, ind_x) = 0;
%         pVisual.pPlot3D = plot3DSolutionAndErrorBehaviour(U_A(2:(end-1), :), Ua(2:(end-1), :), pVisual.pPlot3D, pTask);
        %%
        er1(k)   = max(max(abs(U_E(1:s1:end, 1:s2:end) - U_A(:, :))));
        er1_e(k) = max(max(abs(U_E(1:s1:end, 1:s2:end))));
        er2(k)   = sqrt(sum(sum(abs(U_E(1:s1:end, 1:s2:end) - U_A(:, :)).^2*h*delta)));
        er2_e(k) = sqrt(sum(sum(abs(U_E(1:s1:end, 1:s2:end)).^2*h*delta)));
        %% solution
        %     er1(k)   = max(max(abs(U_E(:, :))));
        %     er1_e(k) = 1;
        %     er2(k)   = sqrt(sum(sum(abs(U_E(:, :)).^2*h*delta)));
        %     er2_e(k) = 1;
        waitbar(k/pCalc.nFreq(end), hw, sprintf('%s %0.2f%s', 'Calculate error: ', k/pCalc.nFreq(end)*100, '%'));
    end
    er1r = er1./er1_e; rError.CNormOfSolution  = er1_e; rError.CAbs  = er1; rError.CRel  = er1r;
    er2r = er2./er2_e; rError.L2NormOfSolution = er2_e; rError.L2Abs = er2; rError.L2Rel = er2r;
    P(pp, 4) = max(er1 );
    P(pp, 5) = max(er2 );
    P(pp, 6) = max(er1r);
    P(pp, 7) = max(er2r);
    close(hw);
    %% Plot absolute and relative errors
    pVisual.pPlotError.sFileSuffix = [ sFileSuffixA '_VS_' sFileSuffixE ];
    if isfield(pTask, 'V_Q')
        pVisual.pPlotError.sFileSuffix = [ 'Q=' num2str(int32(max(pTask.V_Q))) '_' pVisual.pPlotError.sFileSuffix ];
    end
    if isfield(pTask, 'V_lambda')
        pVisual.pPlotError.sFileSuffix = [ 'LAMBDA=' num2str(pTask.V_lambda) '_' pVisual.pPlotError.sFileSuffix ];
    end
    pVisual.pPlotError.sLabelY   = ''; % [ sLabelYE ' vs ' sLabelYA ]; %
    pVisual.pPlotError.nFontSize = pVisual.pPlotError.nFontSize - 2;
    pVisual.pPlotError.nPosition = [0 0 1 1];
    hf = plotErrorBehaviour(pVisual.pPlotError, [ rError.CAbs(pCalc.nFreq); rError.L2Abs(pCalc.nFreq) ]...
        , [ rError.CRel(pCalc.nFreq); rError.L2Rel(pCalc.nFreq) ]); % suptitle([ sLabelYE ' vs ' sLabelYA ]);
    saveImage(hf, pVisual.pPlotError, '', mfilename); close(hf);
    % pVisual.pPlotError.sTitle = {''};
    % hf = plotErrorBehaviour(pVisual.pPlotError, [ rError.CNormOfSolution(pCalc.nFreq); rError.L2NormOfSolution(pCalc.nFreq) ], []);
    % saveImage(hf, pVisual.pPlotError, 'SOLUTION', mfilename); close(hf);
    P(pp, 9) = toc;
end
csvwrite([sExampleName '_' pCalc.method '_' pCalc.type '_' datestr(now,'yyyymmdd_HHMMSS') '.csv'], P);