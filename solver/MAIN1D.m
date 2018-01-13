time = containers.Map({'solve', 'kernel', 'matrix', 'rhs', 'lhs', 'conv', 'error', 'plot', 'file'}, {0, [0 0], 0, 0, 0, 0, 0, 0, 0});
if ( ~pFile.pSolution.bExist || ( pFile.pSolution.bExist && ( ( pCalc.bCalcError && ~pFile.pError.bExist ) ...
        || ( pCalc.bCalcErrorModulus && ~pFile.pErrorModulus.bExist ))))
    %% solve
    tic
    if strcmp(pCalc.type, 'SSP')
        E = calcE(pTask, pCalc);
        pTask.sExampleName = [ '_' pTask.sExampleName ];
    end
    init = psi0(pTask, pInitWave, pCalc.x_mod); % Calc initial function values
    if pCalc.bStoreAllTimeLevels
        U = zeros(pCalc.n_mod, pCalc.m);
        U(1:pCalc.n_mod, 1) = init; % U(pCalc.n_mod, 1)  = 0;
    else
        U = zeros(pCalc.n_mod, 1);
        U(1:pCalc.n_mod) = init; % U(pCalc.n_mod) = 0;
        U_L   = zeros(pCalc.m, 1); U_R    = U_L;
        U_L(1)= U(1);              U_R(1) = U(pCalc.n_mod);
    end
    time('solve') = time('solve') + toc;
    
    %% error
    tic
    switch pCalc.name
        case 'FFDS'
            nn = pCalc.n_perInner; xx = pCalc.xInner;
        case 'FEM'
            nn = pCalc.n_per;      xx = pCalc.x_mod;
    end
    if pCalc.bCalcError || pVisual.bPlotTimeBehaviour % Calc errors
        er1 = zeros(1, pCalc.m); er1_e = er1; er2 = er1; er2_e = er2;
        if pCalc.bStoreAllTimeLevels
            A = U(nn, 1);
        else
            A = U(nn);
        end
        E = transpose(GaussianWave(pTask, pInitWave, pCalc.tau * ( 1 - 1 ), xx));
        outError = calcErrors(A, E, pCalc, pCalc.isConstStep);
        er1(1) = outError.CAbs;  er1_e(1) = outError.CNormOfSolution;
        er2(1) = outError.L2Abs; er2_e(1) = outError.L2NormOfSolution;
    else
        er1 = []; er1_e = er1; er2 = er1; er2_e = er2;
    end
    if pCalc.bCalcErrorModulus
        er1a = zeros(1,pCalc.m); er1a_e = er1a; er2a = er1a; er2a_e = er2a;
        if pCalc.bStoreAllTimeLevels
            UaAbs = abs(U(nn, 1));
        else
            UaAbs = abs(U(nn));
        end
        UeAbs = abs(transpose(GaussianWave(pTask, pInitWave, pCalc.tau * ( 1 - 1 ), xx)));
        outError = calcErrors(UaAbs, UeAbs, pCalc, pCalc.isConstStep);
        er1a(1) = outError.CAbs;  er1a_e(1) = outError.CNormOfSolution;
        er2a(1) = outError.L2Abs; er2a_e(1) = outError.L2NormOfSolution;
    else
        er1a = []; er1a_e = er1a; er2a = er1a; er2a_e = er2a;
    end
    time('error') = time('error') + toc;
    
    %% plot
    tic
    if pVisual.bPlotTimeBehaviour
        if ~pCalc.bCalcError
            E = transpose(GaussianWave(pTask, pInitWave, pCalc.tau * ( 1 - 1 ), xx));
        end
        pVisual.pPlot2D = initGraphics('2D', pVisual.pPlot2D);
        pVisual.pPlot2D = plot2DSolutionBehaviour(A, pVisual.pPlot2D, pTask, 0); % A - E
        saveImage(pVisual.pPlot2D.hFigure, pVisual.pPlot2D, ['k=' num2str(1 - 1)], mfilename); % Save plot of time-level solution
    else
        pVisual.pFrame = [];
    end
    time('plot') = time('plot') + toc;
    
    %% kernel
    tic
    if pCalc.pLBoundary.bUseDTBC % Calc left DTBC's convolution kernel
        pLTask = pTask; pLTask.V_inf = pTask.V_Linf;
        pLCalc = pCalc; pLCalc.name = pCalc.pLBoundary.sName;
        switch pLCalc.name
            case 'FFDS'
                pLCalc.theta = pCalc.pLBoundary.theta;
            case 'FEM'
                pLCalc.N = pCalc.pLBoundary.N;
        end
        [R, c_0] = calcConvolutionKernel(pLTask, pLCalc);
        clear pLTask pLCalc;
    else
        R = []; c_0 =[];
    end
    pCalc.pLBoundary.R = R;
    pCalc.pLBoundary.c_0 = c_0;
    t1 = toc;
    tic
    if pCalc.pRBoundary.bUseDTBC % Calc right DTBC's convolution kernel
        pRTask = pTask; pRTask.V_inf = pTask.V_Rinf;
        pRCalc = pCalc; pRCalc.name = pCalc.pRBoundary.sName;
        switch pRCalc.name
            case 'FFDS'
                pRCalc.theta = pCalc.pRBoundary.theta;
            case 'FEM'
                pRCalc.N = pCalc.pRBoundary.N;
        end
        [R, c_0] = calcConvolutionKernel(pRTask, pRCalc);
        clear pRTask pRCalc;
    else
        R = []; c_0 =[];
    end
    pCalc.pRBoundary.R = R;
    pCalc.pRBoundary.c_0 = c_0;
    t2 = toc;
    time('kernel') = [t1 t2];
    
    %% matrix
    tic
    [pCalc.A, pCalc.Ar] = calcMatrix(pTask, pCalc);
    switch pCalc.method
        case 'TRISYS'
            pCalc.A = full(pCalc.A);
        case 'LU' % get LU decomposition of A
            [ pCalc.LA, pCalc.A ] = lu(pCalc.A); pCalc.LA = pCalc.LA';
        case 'QR' % get QR decomposition of A
            [ pCalc.LA, pCalc.A ] = qr(pCalc.A); pCalc.LA = pCalc.LA';
    end
    time('matrix') = time('matrix') + toc;
    
    %% solve
    tic
    if strcmp(pCalc.type, 'SSP')
        if pCalc.bStoreAllTimeLevels
            V = zeros(pCalc.n_mod, pCalc.m-1); V_L = zeros(pCalc.m-1); V_R = V_L;
        else
            V = zeros(pCalc.n_mod, 1);
        end
    end
    time('solve') = time('solve') + toc;
    
    %%
    for k=2:pCalc.m
        if mod(k-1, 1000) == 0, fprintf('%s Calculate solution on m=%d level...\n', datestr(now, 'HH:MM:SS'), k-1); end
        
        %% solve
        % tSolve = tic;
        ctime = containers.Map({'lhs', 'conv', 'rhs'}, {0, 0, 0});
        if pCalc.bStoreAllTimeLevels
            if strcmp(pCalc.type, 'SSP')
                V(1:pCalc.n_mod, k-1) = E(:) .* U(1:pCalc.n_mod, k-1);
                U(1:pCalc.n_mod, k  ) = E(:) .* calcTimeLevelSolution(V(1:pCalc.n_mod, k-1), transpose(V(1, 1:k-2)), transpose(V(pCalc.n_mod, 1:k-2)), pTask, pCalc, ctime, k);
            else
                [U(1:pCalc.n_mod, k), ctime] = calcTimeLevelSolution(U(1:pCalc.n_mod, k-1), transpose(U(1, 1:k-2)), transpose(U(pCalc.n_mod, 1:k-2)), pTask, pCalc, ctime, k);
                time('lhs')  = time('lhs')  + ctime('lhs');
                time('conv') = time('conv') + ctime('conv');
                time('rhs')  = time('rhs')  + ctime('rhs');
            end
        else
            if strcmp(pCalc.type, 'SSP') % TODO: test it
                V = E(:) .* U; V_L(k-1) = V(1); V_R(k-1) = V(pCalc.n_mod);
                U = E(:) .* calcTimeLevelSolution(V, transpose(V_L(1:k-2)), transpose(V_R(1:k-2)), pTask, pCalc, ctime, k);
            else
                [U, ctime] = calcTimeLevelSolution(U, U_L, U_R, pTask, pCalc, ctime, k);
                time('lhs')  = time('lhs')  + ctime('lhs');
                time('conv') = time('conv') + ctime('conv');
                time('rhs')  = time('rhs')  + ctime('rhs');
                U_L(k) = U(1); U_R(k) = U(pCalc.n_mod);
            end
        end
        % time('solve') = time('solve') + toc(tSolve);
        
        %% error
        tic
        if pCalc.bCalcError || pVisual.bPlotTimeBehaviour
            if pCalc.bStoreAllTimeLevels
                A = U(nn, k);
            else
                A = U(nn);
            end
            E = transpose(GaussianWave(pTask, pInitWave, pCalc.tau * ( k - 1 ), xx));
            outError = calcErrors(A, E, pCalc, pCalc.isConstStep);
            er1(k) = outError.CAbs;  er1_e(k) = outError.CNormOfSolution;
            er2(k) = outError.L2Abs; er2_e(k) = outError.L2NormOfSolution;
        end
        if pCalc.bCalcErrorModulus
            if pCalc.bStoreAllTimeLevels
                UaAbs = abs(U(nn, k));
            else
                UaAbs = abs(U(nn));
            end
            UeAbs = abs(transpose(GaussianWave(pTask, pInitWave, pCalc.tau * ( k - 1 ), xx)));
            outError = calcErrors(UaAbs, UeAbs, pCalc, pCalc.isConstStep);
            er1a(k) = outError.CAbs;  er1a_e(k) = outError.CNormOfSolution;
            er2a(k) = outError.L2Abs; er2a_e(k) = outError.L2NormOfSolution;
        end
        time('error') = time('error') + toc;
        
        %% plot
        tic
        if sum(k-1==linspace(0,pCalc.m-1,pVisual.nPlotCount+1))>0 % mod(k-1,( pCalc.m - 1 ) / pVisual.nPlotCount) == 0 % Log or/and plot the solution
            if pVisual.bPlotTimeBehaviour
                if ~pCalc.bCalcError
                    E = transpose(GaussianWave(pTask, pInitWave, pCalc.tau * ( k - 1 ), xx));
                end
                pVisual.pPlot2D = plot2DSolutionBehaviour(A, pVisual.pPlot2D, pTask, k-1); % A - E
                saveImage(pVisual.pPlot2D.hFigure, pVisual.pPlot2D, ['k=' num2str(k-1)], mfilename);
            end
        end
        time('plot') = time('plot') + toc;
    end
    
    %% file
    tic
    saveVariable(U, pFile.pSolution.sFullName, pFile.pSolution.bSave, mfilename); % Save solution
    time('file') = time('file') + toc;
    
    %% error
    tic
    if pCalc.bCalcError % Save error
        er1r = er1./er1_e;              er2r = er2./er2_e;
        rError.CAbs = er1 ;             rError.L2Abs = er2 ;
        rError.CRel = er1r;             rError.L2Rel = er2r;
        rError.CNormOfSolution = er1_e; rError.L2NormOfSolution = er2_e;
        saveVariable(rError, pFile.pError.sFullName, pFile.pError.bSave, mfilename);
    else
        rError = [];
    end
    if pCalc.bCalcErrorModulus % Save error (modulus of solution)
        er1ar = er1a./er1a_e;           er2ar = er2a./er2a_e;
        rErrorModulus.CAbs = er1a ;     rErrorModulus.L2Abs = er2a;
        rErrorModulus.CRel = er1ar;     rErrorModulus.L2Rel = er2ar;
        saveVariable(rErrorModulus, pFile.pErrorModulus.sFullName, pFile.pErrorModulus.bSave, mfilename);
    else
        rErrorModulus = [];
    end
    time('error') = time('error') + toc;
else
    fprintf('%s Load solution...\n', datestr(now, 'HH:MM:SS'));
    
    %% file
    tic
    U = loadVariable(pFile.pSolution.sFullName, true, mfilename);
    time('file') = time('file') + toc;
    
    %% plot
    tic
    if pVisual.bPlotTimeBehaviour % Plot solution (or error) behaviour
        pVisual.pPlot2D = initGraphics('2D', pVisual.pPlot2D);
        for k=1:pCalc.m
            if sum(k-1==linspace(0,pCalc.m-1,pVisual.nPlotCount+1))>0 % mod(k-1, ( pCalc.m - 1 ) / pVisual.nPlotCount) == 0
                switch pCalc.name
                    case 'FFDS'
                        A = U(pCalc.n_perInner, k);
                    case 'FEM'
                        A = U(pCalc.n_per, k);
                end
                pVisual.pPlot2D = plot2DSolutionBehaviour(A, pVisual.pPlot2D, pTask, k-1);
                saveImage(pVisual.pPlot2D.hFigure, pVisual.pPlot2D, [ 'k=' num2str(k-1) ], mfilename);
            end
        end
    end
    time('plot') = time('plot') + toc;
    
    fprintf('%s Load error...\n', datestr(now, 'HH:MM:SS'));
    %% error
    tic
    if pCalc.bCalcError
        rError = loadVariable(pFile.pError.sFullName, true, mfilename); % Load error
        er1 = rError.CAbs ; er1r = rError.CRel ;
        er2 = rError.L2Abs; er2r = rError.L2Rel;
    else
        rError = [];
    end
    rErrorModulus = [];
    time('error') = time('error') + toc;
end
pCalc.U = U;
time('solve') = time('solve') + (time('lhs')+time('conv')+time('rhs')); time('others') = 0;
time('total') = time('solve') + sum(time('kernel')) + time('matrix');
% time('others') = time('solve') - (time('lhs')+time('conv')+time('rhs'));
T = datestr(datenum(0,0,0,0,0,[ time('total') time('solve') time('lhs') time('conv') time('rhs') ...
    time('others') time('kernel') time('matrix') time('error') time('plot') time('file') ]), 'HH:MM:SS.FFF');
fprintf('%s Timing:\n total=%s\n\t solve =%s\n\t\t lhs   =%s\n\t\t conv  =%s\n\t\t rhs   =%s\n\t\t others=%s\n\t kernel=(%s,%s)\n\t matrix=%s\n error=%s\n plot =%s\n file =%s\n', ...
    datestr(now, 'HH:MM:SS'), T(1, :), T(2, :), T(3, :), T(4, :), T(5, :), T(6, :), T(7, :), T(8, :), T(9, :), T(10, :), T(11, :), T(12, :));
% time('solve'), time('lhs'), time('conv'), time('rhs'), time('solve')-(time('lhs')+time('conv')+time('rhs')), ...
% time('kernel'), time('matrix'), time('error'), time('plot'), time('file'));