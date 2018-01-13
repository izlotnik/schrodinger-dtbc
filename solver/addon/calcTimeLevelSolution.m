function [U, time] = calcTimeLevelSolution(U, U_L, U_R, pTask, pCalc, time, k)
hbar  = pTask.hbar;
B_inf = pTask.B_inf;

%% solve.rhs
tic
b = pCalc.Ar * U;
time('rhs') = time('rhs') + toc;

%% solve.conv, solve.rhs
if pCalc.pLBoundary.bUseDTBC
    tic
    if k > 3
        Convolution = conv2(pCalc.pLBoundary.R(3:k), U_L(1:k-2), 'valid');
    else
        Convolution = 0;
    end
    time('conv') = time('conv') + toc;
    tic
    b(1) = b(1) + hbar^2*B_inf*pCalc.pLBoundary.c_0*Convolution;
    time('rhs') = time('rhs') + toc;
end
if pCalc.pRBoundary.bUseDTBC
    tic
    if k > 3
        Convolution = conv2(pCalc.pRBoundary.R(3:k), U_R(1:k-2), 'valid');
    else
        Convolution = 0;
    end
    time('conv') = time('conv') + toc;
    tic
    n = pCalc.n_mod;
    b(n) = b(n) + hbar^2*B_inf*pCalc.pRBoundary.c_0*Convolution;
    time('rhs') = time('rhs') + toc;
end

%% solve.lhs
tic
switch pCalc.method
    case 'TRISYS'
        U = TRISYS(diag(pCalc.A,-1), diag(pCalc.A), diag(pCalc.A, 1), b);
    case {'LU', 'QR'}
        U = pCalc.A\(pCalc.LA*b);
end
time('lhs') = time('lhs') + toc;