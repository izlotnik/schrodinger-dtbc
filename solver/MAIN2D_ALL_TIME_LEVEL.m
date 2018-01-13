%% Solve task
bCalcError = pCalc.bCalcError; % Backup default value
pCalc.bStoreAllTimeLevels = true;
rCalc.Solution = [];
if ~pFile.pSolution.bExist
    sTask = pTask;
    sTask.B_inf = pTask.B1_inf;
    sTask.sDimension = '1D';
    sTask.nDimension = 1;
    pCalc.bCalcError = false; % Redefine default value
    initFuncVal = calcDST( psi0( pCalc.h_mod*((1:pCalc.n_mod)-1), pCalc.delta*((2:(pCalc.K-1))-1), pInitWave, pTask ) );
    if ~pCalc.pLBoundary.bUseDTBC
        leftFuncVal  = calcDST( GaussianWave2D( 0, pCalc.delta*((2:(pCalc.K-1))-1), pCalc.tau*((1:pCalc.m)-1), pInitWave, pTask ) );
    end
    if ~pCalc.pRBoundary.bUseDTBC
        rightFuncVal = calcDST( GaussianWave2D( pCalc.X_0, pCalc.delta*((2:(pCalc.K-1))-1), pCalc.tau*((1:pCalc.m)-1), pInitWave, pTask ) );
    end
    ll_max = pCalc.K-2;
    U = zeros(ll_max,pCalc.n_mod,pCalc.m); U1 = zeros(pCalc.n_mod,pCalc.m);
    h = waitbar(0, 'Please wait...'); % create waitbar
    for ll=1:ll_max
        lambda = ( ( 2 / pCalc.delta ) * sin ( pi * pCalc.delta * ll / ( 2 * pTask.Y ) ) ) ^ 2;
        mu     = 1 - pCalc.eta * ( pCalc.delta ^ 2 ) * lambda;
        if pCalc.pLBoundary.bUseDTBC
            sTask.V_Linf = pTask.V_Linf + (pTask.hbar^2/2)*pTask.B2_inf*lambda/mu;
            pCalc.pLBoundary.nMultiplier = 1;
        else
            pCalc.leftFuncVal = leftFuncVal(ll,:);
        end
        if pCalc.pRBoundary.bUseDTBC
            sTask.V_Rinf = pTask.V_Rinf + (pTask.hbar^2/2)*pTask.B2_inf*lambda/mu;
            pCalc.pRBoundary.nMultiplier = 1;
        else
            pCalc.rightFuncVal = rightFuncVal(ll,:);
        end
        pCalc.initFuncVal = initFuncVal(ll,:);
        pLTask = sTask; pLCalc = pCalc;
        if pCalc.pLBoundary.bUseDTBC % Calc left DTBC's convolution kernel
            pLTask.V_inf = sTask.V_Linf;
            pLCalc.name = pCalc.pLBoundary.sName;
            switch pLCalc.name
                case 'FFDS'
                    pLCalc.theta = pCalc.pLBoundary.theta;
                case 'FEM'
                    pLCalc.N = pCalc.pLBoundary.N;
            end
            [R, c_0] = calcConvolutionKernel( pLTask, pLCalc );
            c_0 = pCalc.pLBoundary.nMultiplier * c_0;
        else
            R = []; c_0 =[];
        end
        clear pLTask pLCalc;
        pCalc.pLBoundary.R   = R;
        pCalc.pLBoundary.c_0 = c_0;
        pRTask = sTask; pRCalc = pCalc;
        if pCalc.pRBoundary.bUseDTBC % Calc right DTBC's convolution kernel
            pRTask.V_inf = sTask.V_Rinf;
            pRCalc.name = pCalc.pRBoundary.sName;
            switch pRCalc.name
                case 'FFDS'
                    pRCalc.theta = pCalc.pRBoundary.theta;
                case 'FEM'
                    pRCalc.N = pCalc.pRBoundary.N;
            end
            [R, c_0] = calcConvolutionKernel( pRTask, pRCalc );
            c_0 = pCalc.pRBoundary.nMultiplier * c_0;
        else
            R = []; c_0 =[];
        end
        clear pRTask pRCalc;
        pCalc.pRBoundary.R   = R;
        pCalc.pRBoundary.c_0 = c_0;
        switch pCalc.name % Calc matrix and coefficients
            case 'FFDS'
                [Va, Vd, Vc, VaR, VdR, VcR, A, Ar] = calcMatrix( sTask, pCalc );
                vRight = zeros(pCalc.n_mod,3); vRight(1,1:2) = [ VdR(1) VcR(1) ];
                for j=2:(pCalc.n_mod-1)
                    vRight(j,1:3) = [ VaR(j-1) VdR(j) VcR(j) ];
                end
                vRight(pCalc.n_mod,1:2) = [ VaR(pCalc.n_mod-1) VdR(pCalc.n_mod) ];
                pCalc.Va = Va; pCalc.Vc = Vc; pCalc.Vd = Vd;
                pCalc.VaR = VaR; pCalc.VcR = VcR; pCalc.VdR = VdR;
                pCalc.vRight = vRight;
                pCalc.A = A; pCalc.Ar = Ar;
            case 'FEM'
                [vLeft, vRight, A] = calcMatrixFEM( sTask, pCalc );
                if rcond(A) < 10^(-8) % ~ cond(A) > 10^8
                    pCalc.pRBoundary.c_0 = - pCalc.pRBoundary.c_0;
                    pCalc.pLBoundary.c_0 = - pCalc.pLBoundary.c_0;
                    [vLeft, vRight, A] = calcMatrixFEM( sTask, pCalc );
                    disp('Unexpected change of onvolution kernel sign!');
                    keyboard
                end
                pCalc.vLeft = vLeft; pCalc.vRight = vRight;
        end
        switch pCalc.method
            case 'LU' % get LU decomposition of A
                [pCalc.LA, pCalc.UA] = lu(sparse(A)); % [LA, UA, PA] = lu(A);
            case 'QR' % get QR decomposition of A
                [pCalc.LA, pCalc.UA] = qr(sparse(A));
        end
        Ulast = zeros(pCalc.n_mod, 1); Ulast(1:pCalc.n_mod) = pCalc.initFuncVal;
        U_L   = zeros(pCalc.m,1);   U_R    = U_L;
        U_L(1)= Ulast(1);           U_R(1) = Ulast(pCalc.n_mod);
        rCalc.Solution(:,1) = Ulast;
        for k=2:pCalc.m
            Ulast = calcTimeLevelSolution( Ulast, U_L, U_R, sTask, pCalc, k );
            U_L(k) = Ulast(1);
            U_R(k) = Ulast(pCalc.n_mod);
            rCalc.Solution(:,k) = Ulast;
        end
        U(ll,:,:) = rCalc.Solution;
        waitbar(ll/(ll_max), h, sprintf('%s %0.2f%s', 'Calculating (', ll*100/(ll_max), '%)')); % update waitbar
    end
    close(h); % close waitbar
    saveVariable( U, pFile.pSolution.sFullName, pFile.pSolution.bSave, mfilename );
else % Load solution
    U = loadVariable( pFile.pSolution.sFullName, pFile.pSolution.bLoad, mfilename ); U = U.pVariable;
end
pCalc.bCalcError = bCalcError; % Restore default value
if pCalc.bCalcError
    er1 = zeros(1, pCalc.m); er1_e = er1; er1r = er1;
    er2 = er1;               er2_e = er2; er2r = er2;
end
if pCalc.bCalcError || pVisual.bPlotTimeBehaviour
    for k=1:pCalc.m
        Um = calcIDST( U(:,:,k) );
        if pCalc.bCalcError % Calc time level error of solution
            pApproximateSolution = Um(:, pCalc.n_perInner);
            pExactSolution = GaussianWave2D( pCalc.xInner, pCalc.delta*((2:(pCalc.K-1))-1), pCalc.tau*(k-1), pInitWave, pTask );
            rError = calcErrors( pApproximateSolution, pExactSolution, pCalc, true );
            er1(k)   = rError.CAbs;            er2(k)   = rError.L2Abs;
            er1_e(k) = rError.CNormOfSolution; er2_e(k) = rError.L2NormOfSolution;
            er1r(k)  = er1(k)/er1_e(k);        er2r(k)  = er2(k)/er2_e(k);
        end
        if ( mod( k, ( pCalc.m - 1 ) / pVisual.nPlotCount ) == 0 || k == 1 ) && pVisual.bPlotTimeBehaviour % Plot solution (or error) behaviour
            pExactSolution = GaussianWave2D( pCalc.xInner, pCalc.delta*((2:(pCalc.K-1))-1), pCalc.tau*(k-1), pInitWave, pTask );
            if k == 1
                pVisual.pPlot3D = initGraphics( '3D', pVisual.pPlot3D );
            else
                Um(abs(Um) < 2*eps) = 0;
                switch pVisual.pPlot3D.sType
                    case 'modulus'
                        pVisual.pPlot3D.sLabelZ = [ '|\Psi(x,y,t_{'  num2str(k) '})|'];
                    case 'real'
                        pVisual.pPlot3D.sLabelZ = [ 'Re\Psi(x,y,t_{' num2str(k) '})' ];
                    case 'imaginary'
                        pVisual.pPlot3D.sLabelZ = [ 'Im\Psi(x,y,t_{' num2str(k) '})' ];
                end
            end
            pVisual.pPlot3D = plot3DSolutionBehaviour( Um, pVisual.pPlot3D, pTask, k );
            % pVisual.pPlot3D = plot3DSolutionBehaviour( Um - pExactSolution, pVisual.pPlot3D, pTask, k );
            % saveImage( gcf, pVisual.pPlot3D, ['M=' num2str(k)], mfilename );
        end
    end
end