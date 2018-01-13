function E = calcE(pTask, pCalc)
switch pTask.sDimension % set values for dimentional parameters
    case '1D'
        x = pCalc.h_mod*((0:pCalc.n_mod)-1/2);
        t = 1i * pCalc.tau / ( 4 * pTask.hbar ) * ( V(pTask, x) - pTask.V_Rinf ) ./ rho(pTask, x);
        t = ( t(1:(end-1)) + t(2:end) ) / 2; % x averaging
    case '2D' % rho \equiv 1
        x = pCalc.h_mod*((0:pCalc.n_mod)-1/2);
        y = pCalc.delta*((1:pCalc.K-1  )-1/2);
        dV = repmat(transpose(V(pTask, x)), 1, length(y)); % dV=Q\chi(x)
        if isfield(pTask, 'V_y')
            ind = y <= pTask.V_y(1) | y >= pTask.V_y(2);
            if isfield(pCalc, 'subtype') && strcmp(pCalc.subtype, 'CONST')
                dV(:, ind) = 0;
            else
                dV(:, ind) = - dV(:, ind); dV(:, ~ind) = 0;
            end
        end
        dV = ( dV(1:(end-1), :) + dV(2:end, :) ) / 2; % x averaging
        dV = ( dV(:, 1:(end-1)) + dV(:, 2:end) ) / 2; % y averaging
        t = 1i * pCalc.tau / ( 4 * pTask.hbar ) * dV;
end
if isfield(pCalc, 'subtype') && strcmp(pCalc.subtype, 'EXP')
    E = exp( - 2 * t );
else
    E = ( 1 - t ) ./ ( 1 + t );
end
E = transpose(E); % E ~ (Y x X)
end