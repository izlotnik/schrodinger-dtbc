%% Solve task
bCalcError        = pCalc.bCalcError;
bCalcErrorModulus = pCalc.bCalcErrorModulus;
pCalc.bStoreAllTimeLevels = true;
sSaveType = 'SaveAllIntoSingleFile'; % 'SaveByTimeLevel'
if ( ~pFile.pSolution.bExist )
    sTask = pTask;
    sTask.B_inf = pTask.B1_inf;    
    sTask.sDimension = '1D'; 
    sTask.nDimension = 1;    
    pCalc.bCalcError = false;
    pCalc.bCalcErrorModulus = false;
    initFuncVal = calcDST( psi0( pCalc.h_mod*((1:pCalc.n_mod)-1)...
                               , pCalc.delta*((2:(pCalc.K-1))-1), pInitWave, pTask ) );
    if ( ~ pCalc.pLBoundary.bUseDTBC )
        leftFuncVal = calcDST( GaussianWave2D( 0, pCalc.delta*((2:(pCalc.K-1))-1)...
                             , pCalc.tau*((1:pCalc.m)-1), pInitWave ) );
        
    end
    if ( ~ pCalc.pRBoundary.bUseDTBC )
        rightFuncVal = calcDST( GaussianWave2D( pCalc.X_0, pCalc.delta*((2:(pCalc.K-1))-1)...
                              , pCalc.tau*((1:pCalc.m)-1), pInitWave ) );
    end
    ll_max = pCalc.K-2;
    switch sSaveType
        case 'SaveAllIntoSingleFile'
            U = zeros(ll_max,pCalc.n_mod,pCalc.m);
    end
    % create waitbar
    h = waitbar(0, 'Please wait...');
    for ll=1:ll_max
        % update waitbar
        waitbar((ll-1)/(ll_max), h, sprintf( '%s %0.2f%s' ...
            , 'Calculating (', (ll-1)*100/(ll_max), '%)' ) );
        lambda = ( ( 2 / pCalc.delta )...
               * sin ( pi * pCalc.delta * ll / ( 2 * pTask.Y ) ) ) ^ 2;
        mu     = 1 - pCalc.eta * ( pCalc.delta ^ 2 ) * lambda;
        if ( pCalc.pLBoundary.bUseDTBC )
            sTask.V_Linf = pTask.V_Linf + (pTask.hbar^2/2)*pTask.B2_inf*lambda/mu;
            pCalc.pLBoundary.nMultiplier = 1;
        else
            pCalc.pLBoundary.nMultiplier = 1;
            pCalc.leftFuncVal = leftFuncVal(ll,:);
        end
        if ( pCalc.pRBoundary.bUseDTBC )
            sTask.V_Rinf = pTask.V_Rinf + (pTask.hbar^2/2)*pTask.B2_inf*lambda/mu;
            pCalc.pRBoundary.nMultiplier = 1;
        else
            pCalc.pRBoundary.nMultiplier = 1;
            pCalc.rightFuncVal = rightFuncVal(ll,:);
        end
        pCalc.initFuncVal = initFuncVal(ll,:);
        [rCalc, pVisual, rError, rModulError] = SOLVER( sTask, pCalc, pVisual, pInitWave );
        switch sSaveType
            case 'SaveAllIntoSingleFile'
                U(ll,:,:) = rCalc.Solution;
            case 'SaveByTimeLevel'
                tic
                for mm=1:pCalc.m
                    dlmwrite( [pCalc.sWorkDirectory '\FOURIER_COEF_M=' num2str(mm) '.csv']...
                            , transpose(rCalc.Solution(:,mm))...
                            , '-append', 'delimiter', '\t', 'precision', 16 );
                end
                toc
        end
        % Plot kernel
        % plotkernel
        % update waitbar
        waitbar(ll/(ll_max),h, sprintf( '%s %0.2f%s' ...
            , 'Calculating (', ll*100/(ll_max), '%)' ) );
    end
    % close waitbar
    close(h);
    switch sSaveType
        case 'SaveAllIntoSingleFile'
            % Save Fourier coefficients of solution
            saveVariable ( U, pFile.pSolution.sFullName, pFile.pSolution.bSave, mfilename );
    end
else
    switch sSaveType
        case 'SaveAllIntoSingleFile'
            % Load Fourier coefficients of solution
            U = loadVariable ( pFile.pSolution.sFullName, pFile.pSolution.bLoad, mfilename );
            U = U.pVariable;
    end
end
pCalc.bCalcError        = bCalcError;
pCalc.bCalcErrorModulus = bCalcErrorModulus;
if ( pCalc.bCalcError )
    er1 = zeros(1, pCalc.m); er1_e = er1; er1r = er1;
    er2 = er1;               er2_e = er2; er2r = er2;
end
if ( pCalc.bCalcErrorModulus )
    er1a = zeros(1, pCalc.m); er1a_e = er1a; er1ar = er1a;
    er2a = er1a;              er2a_e = er2a; er2ar = er2a;
end
for mm=1:pCalc.m
    % Calc inverse discrete sine transform of Fourier coefficients
    switch sSaveType
        case 'SaveAllIntoSingleFile'
            Um = calcIDST( U(:,:,mm) );
        case 'SaveByTimeLevel'
            Um = calcIDST(dlmread( [pCalc.sWorkDirectory '\FOURIER_COEF_M=' num2str(mm) '.csv']...
                                 ,  '\t'));
    end
    if ( pCalc.bCalcError )
        pApproximateSolution = Um(:, pCalc.n_perInner);
        pExactSolution = GaussianWave2D( pCalc.xInner...
                         , pCalc.delta*((2:(pCalc.K-1))-1), pCalc.tau*(mm-1), pInitWave );
        % Calc time level error of solution
        rError = calcErrors( pApproximateSolution, pExactSolution, pCalc, true );
        er1(mm)   = rError.CAbs;
        er1_e(mm) = rError.CNormOfSolution;
        er2(mm)   = rError.L2Abs;
        er2_e(mm) = rError.L2NormOfSolution;
        er1r(mm)  = er1(mm)/er1_e(mm);
        er2r(mm)  = er2(mm)/er2_e(mm);
    end
    if ( pCalc.bCalcErrorModulus )
        pApproximateSolution = Um(:, pCalc.n_perInner);
        pExactSolution = GaussianWave2D( pCalc.xInner...
                         , pCalc.delta*((2:(pCalc.K-1))-1), pCalc.tau*(mm-1), pInitWave );
        % Calc time level error of modulus of solution
        rError = calcErrors( abs(pApproximateSolution), abs(pExactSolution), pCalc, true );
        er1a(mm)   = rError.CAbs;
        er1a_e(mm) = rError.CNormOfSolution;
        er2a(mm)   = rError.L2Abs;
        er2a_e(mm) = rError.L2NormOfSolution;
        er1ar(mm)  = er1a(mm)/er1a_e(mm);
        er2ar(mm)  = er2a(mm)/er2a_e(mm);
    end
   if ( ( mod( mm, ( pCalc.m - 1 ) / pVisual.nPlotCount ) == 0 || mm == 1 )...
           && pVisual.bPlotTimeBehaviour )
        switch pVisual.pPlot3D.sType
            case 'modulus'
                pVisual.pPlot3D.sLabelZ = [ '|\Psi(x,y,t_{'     num2str(mm) '})|' ];
            case 'real'
                pVisual.pPlot3D.sLabelZ = [ 'Re\Psi(x,y,t_{'    num2str(mm) '})'  ];
            case 'imaginary'
                pVisual.pPlot3D.sLabelZ = [ 'Im\Psi(x,y,t_{'    num2str(mm) '})'  ];
        end
        pExactSolution = GaussianWave2D( pCalc.h_mod*((1:pCalc.n_mod)-1)...
                         , pCalc.delta*((2:(pCalc.K-1))-1), pCalc.tau*(mm-1), pInitWave );
        Um( abs(Um) < 2*eps ) = 0;
        if ( mm == 1 )
            pVisual.pPlot3D = initGraphics( '3D', pVisual.pPlot3D );
            pVisual.pPlot3D = plot3DSolutionBehaviour( Um, pVisual.pPlot3D, pTask, mm );
            % plot error
            % pVisual.pPlot3D = plot3DSolutionBehaviour( Um - pExactSolution, pVisual.pPlot3D, pTask );
        else
            pVisual.pPlot3D = plot3DSolutionBehaviour( Um, pVisual.pPlot3D, pTask, mm );
            % plot error
            % pVisual.pPlot3D = plot3DSolutionBehaviour( Um - pExactSolution, pVisual.pPlot3D, pTask );
        end
        % Save image
        saveImage ( gcf, pVisual.pPlot3D, ['M=' num2str(mm)], mfilename );
    end
end