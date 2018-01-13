clearContext; format short; echo on;
%% Set base parameters
sExampleName   = 'EX22a';
sDataDirectory = 'C:\MATLAB\R2011b'; % 'C:\MATLAB'; % 'C:\MATLAB\R2009b\work'; % 'H:\MATLAB-DATA'; %
sBaseRunID     = '';             % 'LAST'; % '20100205_170006_125'; %
%% Set task parameters
[pTask pInitWave] = setTaskParams('Schrodinger', '2D', sExampleName, sDataDirectory, sBaseRunID); pTaskE = pTask;
echo off; tic
P = [ ...
     2*1200 2*600 4*64 ... % EX22a
%    1*700 1*800 1*128 ... % EX02m
%    1*(7/12)*1200/2 2*(2/3)*600/1 4*64/1 
%    1*1200/1 1*(3/5)*1000/1 1*64/1
%    1*1200/4 1*(3/5)*1000/3 1*128/4 
%    1*1200/1 1*(3/5)*1000/1 1*128/2
%    1*1200/2 1*(3/5)*1000/2 1*128/4
%    1*1200/4/10 1*(3/5)*1000/2 1*128/4
    ];
% s = size(P, 1)+1;
% for j=2.^(-2:1)
%     for m=2.^(-2:1)
%         for k=2.^(-2:1)
%             % EX02m         
%             P(s, :) = [ j*(7/12)*1200 m*800 k*128 ];
%             s = s + 1;
%         end
%     end
% end
echo off; tic
if true % false % false - plot solution ONLY
    for pp=1:size(P, 1)
        pTask = pTaskE; pCalc = []; pFile = []; pVisual = []; pp
        %% Set calculation parameters
        [pCalc   pFile] = setCalcParams(pTask, {'FFDS'}, [1/12 P(pp, 1) P(pp, 2) NaN 0 1/12 P(pp, 3)], {'TRISYS'}); % , 'SSP', 'CONST'
        %% Set visualization parameters
        [pVisual pCalc] = setVisualParams(pTask, pCalc, [1 1 50 0 0], [0 0 0 0], [0 1]); % plot sloution and error % [1 10 0 0]
        %% Solve task
        % pVisual.bPlotTimeBehaviour = false;
        % pCalc.bCalcError             = false; % calc error
        % pVisual.pPlot2D.bSaveFigure  = false;
        % pVisual.pPlot3D.bPlayPause   = true ; % view every figure via KEYBOARD
        MAIN2D_ONE_TIME_LEVEL % MAIN2D_ALL_TIME_LEVEL % MAIN2D %
        P(pp, 4) = max(er1 );
        P(pp, 5) = max(er2 );
        P(pp, 6) = max(er1r);
        P(pp, 7) = max(er2r);
        P(pp, 8) = toc;
%         for k=pCalc.nFreq
%             if k > 1
%                 sFullName = [ pCalc.sWorkDirectory '\' pFile.pSolution.sName '_M=' num2str((k-1)) '.' pFile.pSolution.sExtension ];
%             else
%                 sFullName = pFile.pSolution.sFullName;
%             end
%             U_A = loadVariable(sFullName, pFile.pSolution.bLoad, mfilename);
%             U_E = GaussianWave(pTask, pInitWave, pCalc.tau*(k-1), pCalc.xInner, pCalc.delta*((2:(pCalc.K-1))-1));
%             % U_E = psi0(pCalc.h_mod*((1:pCalc.n_mod)-1), pCalc.delta*((2:(pCalc.K-1))-1), pInitWave, pTask);
%             if k==1
%                 pVisual.pPlot3D = initGraphics('3D', pVisual.pPlot3D);
%             end
%             pVisual.pPlot3D = plot3DSolutionAndErrorBehaviour(U_A, U_A-U_E, pVisual.pPlot3D, pTask);
%         end
        % close(gcf);
        if ~isempty(rError) % Print to command line max error values
            printMaxErrors(rError, [], pCalc);
        end
        if pVisual.bPlotErrorTimeBehaviour && ~isempty(rError) % Plot and save plot time-level error behaviour
            plotErrorBehaviour(pVisual.pPlotError, [ rError.CAbs(pCalc.nFreq); rError.L2Abs(pCalc.nFreq) ], [ rError.CRel(pCalc.nFreq); rError.L2Rel(pCalc.nFreq) ]);
            saveImage(gcf, pVisual.pPlotError, '', mfilename); % close(gcf);
        end
    end
    csvwrite([sExampleName '_' pCalc.method '_' pCalc.type '_' datestr(now,'yyyymmdd_HHMMSS') '.csv'], P);
else
    if pVisual.bPlotTimeBehaviour % Plot solution behaviour
        pVisual.pPlot3D = initGraphics('3D', pVisual.pPlot3D);
        for k=1:pCalc.m
            if mod(k-1, (pCalc.m - 1) / pVisual.pPlot3D.nPlotCount) == 0 && pVisual.bPlotTimeBehaviour % Plot solution behaviour
                Ue = GaussianWave2D(pTask, pInitWave, pCalc.tau*(k-1), pCalc.xInner, pCalc.delta*((2:(pCalc.K-1))-1));
                pVisual.pPlot3D = plot3DSolutionBehaviour(Ue, pVisual.pPlot3D, pTask, k);
                saveImage(gcf, pVisual.pPlot3D, ['M=' num2str(k-1)], mfilename);
            end
        end
    end
end
toc; echo on;
if true
    max(max(abs(GaussianWave(pTask, pInitWave, pCalc.tau*((1:pCalc.m)-1), 0,            pCalc.delta*(1:pCalc.K-2)   )))) % eps0
    max(    abs(GaussianWave(pTask, pInitWave, 0,                         pCalc.X_0,    pCalc.delta*(1:pCalc.K-2)   )))
    max(max(abs(GaussianWave(pTask, pInitWave, pCalc.tau*(pCalc.m-1),     pCalc.xInner, pCalc.delta*([0 pCalc.K-1]) )))) % eps1
    max(max(abs(GaussianWave(pTask, pInitWave, pCalc.Tmax,                pCalc.xInner, pCalc.delta*(1:pCalc.K-2)   )))) % max last value
end