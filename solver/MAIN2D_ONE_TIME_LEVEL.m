%% Solve task
if ~pFile.pSolution.bExist % || ( pFile.pSolution.bExist && ( ( pCalc.bCalcError && ~pFile.pError.bExist ) ) )
    initFuncVal = calcDST(psi0(pTask, pInitWave, pCalc.h_mod*((1:pCalc.n_mod)-1), pCalc.delta*((2:(pCalc.K-1))-1))); % Calc initial function values
    pTask.B_inf = pTask.B1_inf;
    sTask = pTask;
    if strcmp(pCalc.type, 'SSP')
       E = calcE(pTask, pCalc);
       % I_x = pTask.V_x/pTask.X*(pCalc.n-1); I_y = pTask.V_y/pTask.Y*(pCalc.K-1); imag(E(I_x(1):(I_x(2)+2),[1:(I_y(1)+1) (I_y(2)-1):end]))
       figure();
       x = pCalc.h_mod*((0:(pCalc.n_mod-1))); y = pCalc.delta*(1:(pCalc.K-2));
       [X, Y] = meshgrid(x,y);
       surfc(X, Y, real(E), 'FaceColor', [0 1 0], 'FaceAlpha', 0.10, 'EdgeAlpha', 1); view(2); rotate3d on;
       xlim([x(1) x(end)]); ylim([y(1) y(end)]);
       if isfield(pCalc, 'subtype') && strcmp(pCalc.subtype, 'CONST')
           sTask.sExampleName = [ '_' sTask.sExampleName ];
       end
    end
    sTask.sDimension = '1D';
    ll_max = pCalc.K-2;
    K_L    = zeros(ll_max, pCalc.m); K_R = K_L;
    c_L    = zeros(ll_max,       1); c_R = c_L;
    switch pCalc.method
        case {'LU','QR'}
            LA = zeros(ll_max, pCalc.n_mod, pCalc.n_mod); UA = LA;
            pCalc.LA = zeros(pCalc.n_mod, pCalc.n_mod); pCalc.UA = pCalc.LA;
    end
    Va     = zeros(ll_max, pCalc.n_mod-1); Vc = Va; Vd  = zeros(ll_max, pCalc.n_mod); VaR = Va; VdR = Vd; VcR = Vc;
    vRight = zeros(ll_max, pCalc.n_mod,3); pCalc.vRight = zeros(pCalc.n_mod, 3); S = zeros(ll_max, 1);
    Ulast  = zeros(pCalc.n_mod, ll_max);   U_L = zeros(ll_max, pCalc.m); U_R = U_L;
    hw = waitbar(0, 'Calculate solution: please wait...');
    for ll=1:ll_max
        lambda = ( ( 2 / pCalc.delta ) * sin( pi * pCalc.delta * ll / ( 2 * pTask.Y ) ) )^2;
        S(ll) = ( pTask.hbar^2 / 2 ) * pTask.B2_inf * lambda / ( 1 - pCalc.eta * pCalc.delta^2 * lambda );
        if pCalc.pLBoundary.bUseDTBC % Calc left DTBC's convolution kernel
            pCalc.pLBoundary.nMultiplier = 1;
            pLTask = sTask; pLTask.V_inf = pTask.V_Linf + S(ll);
            pLCalc = pCalc; pLCalc.name  = pCalc.pLBoundary.sName;
            switch pLCalc.name
                case 'FFDS'
                    pLCalc.theta = pCalc.pLBoundary.theta;
                case 'FEM'
                    pLCalc.N     = pCalc.pLBoundary.N;
            end
            [ K_L(ll,:), c_L(ll) ] = calcConvolutionKernel(pLTask, pLCalc);
            c_L(ll) = pCalc.pLBoundary.nMultiplier * c_L(ll);
            pCalc.pLBoundary.R   = K_L(ll, :);
            pCalc.pLBoundary.c_0 = c_L(ll);
            clear pLTask pLCalc;
        end
        if pCalc.pRBoundary.bUseDTBC % Calc right DTBC's convolution kernel
            pCalc.pRBoundary.nMultiplier = 1;
            pRTask = sTask; pRTask.V_inf = pTask.V_Rinf + S(ll);
            pRCalc = pCalc; pRCalc.name  = pCalc.pRBoundary.sName;
            switch pRCalc.name
                case 'FFDS'
                    pRCalc.theta = pCalc.pRBoundary.theta;
                case 'FEM'
                    pRCalc.N     = pCalc.pRBoundary.N;
            end
            [ K_R(ll,:), c_R(ll) ] = calcConvolutionKernel(pRTask, pRCalc);
            c_R(ll) = pCalc.pRBoundary.nMultiplier * c_R(ll);
            pCalc.pRBoundary.R   = K_R(ll, :);
            pCalc.pRBoundary.c_0 = c_R(ll);
            clear pRTask pRCalc;
        end
        %% ONLY FOR STEP-WISED POTENTIAL
        sTask.V_Linf = pTask.V_Linf + S(ll);
        if isfield(pTask, 'V_Q')
            sTask.V_Q = pTask.V_Q   + S(ll);
        end
        sTask.V_Rinf = pTask.V_Rinf + S(ll);
        switch pCalc.name % Calc matrix and coefficients
            case 'FFDS'
                [ Va(ll, :), Vd(ll, :), Vc(ll, :), VaR(ll, :), VdR(ll, :), VcR(ll, :), A, ~ ] = calcMatrix(sTask, pCalc);
                vRight(ll, 1, 1:2) = [ VdR(ll, 1) VcR(ll, 1) ];
                for j=2:(pCalc.n_mod-1), vRight(ll, j, 1:3) = [ VaR(ll, j-1) VdR(ll, j) VcR(ll, j) ]; end
                vRight(ll, pCalc.n_mod, 1:2) = [ VaR(ll, pCalc.n_mod-1) VdR(ll, pCalc.n_mod) ];
            case 'FEM'
                [ vLeft, vRight(ll,:,:), A ] = calcMatrixFEM(sTask, pCalc);
                if rcond(A) < 10^(-8) % ~ cond(A) > 10^8
                    pCalc.pRBoundary.c_0 = - pCalc.pRBoundary.c_0;
                    pCalc.pLBoundary.c_0 = - pCalc.pLBoundary.c_0;
                    [ vLeft, vRight(ll,:,:), A ] = calcMatrixFEM(sTask, pCalc);
                    disp('Unexpected sign change of the convolution kernel!'); keyboard
                end
        end
        switch pCalc.method
            case 'LU' % get LU decomposition of A
                [ LA(ll,:,:), UA(ll,:,:) ] = lu(sparse(A)); % [LA, UA, PA] = lu(A);
            case 'QR' % get QR decomposition of A
                [ LA(ll,:,:), UA(ll,:,:) ] = qr(sparse(A));
        end
        Ulast(1:pCalc.n_mod, ll) = initFuncVal(ll, :);
        U_L(ll, 1) = Ulast(          1, ll);
        U_R(ll, 1) = Ulast(pCalc.n_mod, ll);
    end
    if pCalc.bCalcError
        er1 = zeros(1, pCalc.m); er1_e = er1; er1r = er1;
        er2 = er1;               er2_e = er2; er2r = er2;
    end
    U = calcIDST(transpose(Ulast));
    saveVariable(U, pFile.pSolution.sFullName, pFile.pSolution.bSave, mfilename);
    Ua = U(:, pCalc.n_perInner);
    if pCalc.bCalcError % Calc time level error of solution
        Ue = GaussianWave(pTask, pInitWave, pCalc.tau*(1-1), pCalc.xInner, pCalc.delta*((2:(pCalc.K-1))-1));
        rError = calcErrors(Ua, Ue, pCalc, true);
        er1(1) = rError.CAbs;  er1_e(1) = rError.CNormOfSolution;  er1r(1) = er1(1)/er1_e(1);
        er2(1) = rError.L2Abs; er2_e(1) = rError.L2NormOfSolution; er2r(1) = er2(1)/er2_e(1);
    end
    if pVisual.bPlotTimeBehaviour % Plot solution (or error) behaviour
        if ~pCalc.bCalcError
            Ue = GaussianWave(pTask, pInitWave, pCalc.tau*(1-1), pCalc.xInner, pCalc.delta*((2:(pCalc.K-1))-1));
        end
        pVisual.pPlot3D = initGraphics('3D', pVisual.pPlot3D); % Ua(abs(Ua) < 2*eps) = min(min(Ua));
        pVisual.pPlot3D = plot3DSolutionBehaviour(Ua, pVisual.pPlot3D, pTask); % Ua - Ue
        saveImage(pVisual.pPlot3D.hFigure, pVisual.pPlot3D, ['M=' num2str(0)], mfilename);
    end
    if strcmp(pCalc.type, 'SSP')
        fUc = zeros(pCalc.n_mod, ll_max);
        fUc_L = zeros(ll_max, pCalc.m-1); fUc_R = fUc_L;
        fUt_L = zeros(ll_max, pCalc.m  ); fUt_R = fUt_L;
        for ll=1:ll_max
            fUt_L(ll, 1) = initFuncVal(ll, 1); fUt_R(ll, 1) = initFuncVal(ll, pCalc.n_mod);
        end
        % fUt_L(ll, 1) = 0; fUt_R(ll, 1) = 0;
    end
    for k=2:pCalc.m
        if strcmp(pCalc.type, 'SSP')
            Uc  = E .* U;
            fUc = calcDST(Uc);
            fUt = zeros(pCalc.n_mod, ll_max);
        end
        for ll=1:ll_max
            pCalc.pLBoundary.R   = K_L(ll, :);
            pCalc.pLBoundary.c_0 = c_L(ll);
            pCalc.pRBoundary.R   = K_R(ll, :);
            pCalc.pRBoundary.c_0 = c_R(ll);
            pCalc.vRight(:)      = vRight(ll, :, :);
            switch pCalc.method
                case {'LU', 'QR'} % TODO: it's VERY slow - speedup performance
                    pCalc.LA(:) = LA(ll, :, :);
                    pCalc.UA(:) = UA(ll, :, :);
                case 'TRISYS'
                    pCalc.Va = Va(ll, :); pCalc.VaR = VaR(ll, :);
                    pCalc.Vd = Vd(ll, :); pCalc.VdR = VdR(ll, :);
                    pCalc.Vc = Vc(ll, :); pCalc.VcR = VcR(ll, :);
            end
            if strcmp(pCalc.type, 'SSP')
                % fUc_L(ll, k-1) = fUc(ll, 1); fUc_R(ll, k-1) = fUc(ll, pCalc.n_mod);
                % fUt(:, ll) = calcTimeLevelSolution(transpose(fUc(ll, :)), fUc_L(ll, 1:k-2), fUc_R(ll, 1:k-2), pTask, pCalc, k);
                fUt(:, ll) = calcTimeLevelSolution(transpose(fUc(ll, :)), fUt_L(ll, 1:k-2), fUt_R(ll, 1:k-2), pTask, pCalc, k);
                fUt_L(ll, k) = fUt(1, ll); fUt_R(ll, k) = fUt(pCalc.n_mod, ll);
            else
                Ulast(:, ll) = calcTimeLevelSolution(Ulast(:, ll), U_L(ll, :), U_R(ll, :), pTask, pCalc, k);
                U_L( ll,  k) = Ulast(          1, ll);
                U_R( ll,  k) = Ulast(pCalc.n_mod, ll);
            end
        end
        if strcmp(pCalc.type, 'SSP')
            Ut = calcIDST(transpose(fUt));
            U  = E .* Ut;
        else
            U = calcIDST(transpose(Ulast));
        end
        if sum(k==pCalc.nFreq)>0
            sFullName = [ pCalc.sWorkDirectory '\' pFile.pSolution.sName '_M=' num2str(k-1) '.' pFile.pSolution.sExtension ];
            saveVariable(U, sFullName, pFile.pSolution.bSave, mfilename);
        end
        if sum(k==pCalc.nFreq)>0 && pCalc.bCalcError % Calc time level error of solution
            Ua = U(:, pCalc.n_perInner);
            Ue = GaussianWave(pTask, pInitWave, pCalc.tau*(k-1), pCalc.xInner, pCalc.delta*((2:(pCalc.K-1))-1));
            rError = calcErrors(Ua, Ue, pCalc, true);
            er1(k) = rError.CAbs;  er1_e(k) = rError.CNormOfSolution;  er1r(k) = er1(k)/er1_e(k);
            er2(k) = rError.L2Abs; er2_e(k) = rError.L2NormOfSolution; er2r(k) = er2(k)/er2_e(k);
        end
        if mod(k-1, ( pCalc.m - 1 ) / pVisual.nPlotCount) == 0 && pVisual.bPlotTimeBehaviour % Plot solution (or error) behaviour
            if ~pCalc.bCalcError
                Ue = GaussianWave(pTask, pInitWave, pCalc.tau*(k-1), pCalc.xInner, pCalc.delta*((2:(pCalc.K-1))-1));
            end
            Ua = U(:, pCalc.n_perInner); % Ua(abs(Ua) < 2*eps) = min(min(Ua));
            switch pVisual.pPlot3D.sType
                case 'modulus'
                    pVisual.pPlot3D.sLabelZ = [ '|\Psi(x,y,t_{'  num2str(k) '})|'];
                case 'real'
                    pVisual.pPlot3D.sLabelZ = [ 'Re\Psi(x,y,t_{' num2str(k) '})' ];
                case 'imaginary'
                    pVisual.pPlot3D.sLabelZ = [ 'Im\Psi(x,y,t_{' num2str(k) '})' ];
            end
            pVisual.pPlot3D = plot3DSolutionBehaviour(Ua, pVisual.pPlot3D, pTask); % Ua - Ue
            saveImage(pVisual.pPlot3D.hFigure, pVisual.pPlot3D, ['M=' num2str(k-1)], mfilename);
        end
        waitbar(k/pCalc.m, hw, sprintf('%s: %0.2f%s', 'Calculate solution', k/pCalc.m*100, '%'));
    end
    close(hw);
    if pCalc.bCalcError % Save error
        rError.CAbs = er1 ; rError.L2Abs = er2 ;
        rError.CRel = er1r; rError.L2Rel = er2r;
        saveVariable(rError, pFile.pError.sFullName, pFile.pError.bSave, mfilename);
    else
        rError = [];
    end
else
    if pVisual.bPlotTimeBehaviour % Plot solution (or error) behaviour
        pVisual.pPlot3D = initGraphics('3D', pVisual.pPlot3D);
        for k=1:pCalc.m
            if  sum(k==pCalc.nFreq)>0 && mod(k-1, ( pCalc.m - 1 ) / pVisual.pPlot3D.nPlotCount) == 0
                if k > 1
                    sFullName = [ pCalc.sWorkDirectory '\' pFile.pSolution.sName '_M=' num2str(k-1) '.' pFile.pSolution.sExtension ];
                else
                    sFullName = pFile.pSolution.sFullName;
                end
                U = loadVariable(sFullName, pFile.pSolution.bLoad, mfilename);
                switch pVisual.pPlot3D.sType
                    case 'modulus'
                        pVisual.pPlot3D.sLabelZ = [ '|\Psi(x,y,t_{'  num2str(k) '})|'];
                    case 'real'
                        pVisual.pPlot3D.sLabelZ = [ 'Re\Psi(x,y,t_{' num2str(k) '})' ];
                    case 'imaginary'
                        pVisual.pPlot3D.sLabelZ = [ 'Im\Psi(x,y,t_{' num2str(k) '})' ];
                end
                % x = pCalc.h_mod*(0:(pCalc.n_mod-1)); ind_x = x >= pTask.V_x(1) & x <= pTask.V_x(2);
                % y = pCalc.delta*(1:(pCalc.K-2    )); ind_y = y >= pTask.V_y(1) & y <= pTask.V_y(2);
                % U(ind_y, :) = 0; U(:, ~ind_x) = 0; max(max(abs(U)))
                pVisual.pPlot3D = plot3DSolutionBehaviour(U, pVisual.pPlot3D, pTask);
                saveImage(pVisual.pPlot3D.hFigure, pVisual.pPlot3D, ['M=' num2str(k-1)], mfilename);
            end
        end
    end
    if pCalc.bCalcError
        rError = loadVariable(pFile.pError.sFullName, pFile.pError.bLoad, mfilename);
        er1 = rError.CAbs ; er1r = rError.CRel ;              
        er2 = rError.L2Abs; er2r = rError.L2Rel;
    end
    rErrorModulus = [];
end