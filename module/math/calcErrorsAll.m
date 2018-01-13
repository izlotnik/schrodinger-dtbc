function rError = calcErrorsAll(U, pTask, pInitWave, pCalc, params)
if params.isKey('etalon') && ( strcmp(params('etalon'), 'global') || strcmp(params('etalon'), 'local') )
    useEtalon = true;
    U_E = params('U');
    pCalc_E = params('params');
else
    useEtalon = false;
end
if params('calc')
    er1 = zeros(1, pCalc.m); er1_e = er1; er2 = er1; er2_e = er2;
    switch pCalc.name
        case 'FFDS'
            nn = pCalc.n_perInner; xx = pCalc.xInner;
        case 'FEM'
            nn = pCalc.n_per;      xx = pCalc.x_mod;
    end
    if useEtalon
        s = (size(U_E, 2)-1)/(pCalc.m-1); if s~=int32(s), keyboard, end
        x = (pCalc.n_per-1)'*pCalc.h_mod;
    end
    for k=pCalc.nFreq % 1:pCalc.m % 
        if useEtalon
            if pCalc_E.n == pCalc.n && pCalc_E.N == 1 && ~isfield(pCalc, 'N')
                y = U_E(:, (k-1)*s+1);
            else
                D = PFD(U_E(:, (k-1)*s+1), pCalc_E.n-1, pCalc_E.N);
                y = PN(D, pCalc_E.N, pCalc_E.h, x, pCalc_E.n-1);
            end
        else
            y = transpose(GaussianWave(pTask, pInitWave, pCalc.tau*(k-1), xx));
        end
        rError = calcErrors(U(nn, k), y, pCalc, pCalc.isConstStep);
        er1(k) = rError.CAbs ; er1_e(k) = rError.CNormOfSolution ;
        er2(k) = rError.L2Abs; er2_e(k) = rError.L2NormOfSolution;
    end
    if useEtalon % Initial function error is zero (or move it to the plotting function - ?)
        er1(1) = 0; er2(1) = 0;
    end
    rError.CAbs  = er1; rError.CNormOfSolution  = er1_e;
    rError.L2Abs = er2; rError.L2NormOfSolution = er2_e;
    rError.CRel  = er1./er1_e;
    rError.L2Rel = er2./er2_e;
else
    rError = [];
end
end