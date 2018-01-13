function pError = calcErrors(U_A, U_E, pCalc, bIsStepConstant)
h = pCalc.h;
if size(U_A, 2) == 1 % calc 1D error
    pError.CAbs = max(abs(U_E - U_A));
    pError.CNormOfSolution = max(abs(U_E));
    if bIsStepConstant
        switch pCalc.name
            case 'FFDS'
                pError.L2Abs = sqrt(sum(abs(U_E - U_A).^2*h));
                pError.L2NormOfSolution = sqrt(sum(abs(U_E).^2*h));            
            case 'FEM'
                N = pCalc.N;
                n = pCalc.n;
                pError.L2Abs = sqrt(calcNewtonCotesQuadrature(abs(U_E - U_A).^2, (n-1)*N, h, N));
                pError.L2NormOfSolution = sqrt(calcNewtonCotesQuadrature(abs(U_E).^2, (n-1)*N, h, N));
        end
    else
        nInner = pCalc.nInner;
        er1   = abs(U_E(1) - U_A(1));
        er1_e = abs(U_E(1));
        er2   = 0;
        er2_e = 0;
        for j = 1:(nInner+1)
            nm = abs(U_E(j) - U_A(j));
            if nm > er1
                er1 = nm;
            end
            nm2 = abs(U_E(j));
            if nm2 > er1_e
                er1_e = nm2;
            end
            er2   = er2   + abs(U_E(j) - U_A(j))^2*h;
            er2_e = er2_e + abs(U_E(j))^2*h;
        end
        er2   = sqrt(er2);
        er2_e = sqrt(er2_e);
        pError.CAbs  = er1; pError.CNormOfSolution  = er1_e;
        pError.L2Abs = er2; pError.L2NormOfSolution = er2_e;
    end
else % calc 2D error
    delta = pCalc.delta;
    pError.CAbs = max(max(abs(U_E - U_A)));
    pError.CNormOfSolution = max(max(abs(U_E)));
    pError.L2Abs = sqrt(sum(sum(abs(U_E - U_A).^2*h*delta)));
    pError.L2NormOfSolution = sqrt(sum(sum(abs(U_E).^2*h*delta)));
end