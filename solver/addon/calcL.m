function L = calcL(m, param1, param2)
if nargin <= 1
    error(message('calcL:TooFewInputs'));
else if nargin == 2
        phi = param1;
        switch length(phi)
            case 1
                kappa = - exp(1i*phi/2);
                mu    =   cos(   phi/2);
            case 2
                kappa = - exp(1i*(phi(1) + phi(2))/2);
                mu    =   cos(   (phi(1) - phi(2))/2);
            otherwise
                error(message('calcL:InvalidData'));
        end
    else if nargin == 3
            kappa = param1;
            mu    = param2;
        else
            error(message('calcL:TooMuchInputs'));
        end
    end
end
L = NaN * ones(m, 1);
L(1) = 1;
L(2) = kappa * mu;
for k=2:(m-1)
    L(k+1) = ( ( 2 * k - 1 ) / k ) * kappa * mu * L(k) - ( ( k - 1 ) / k ) * kappa^2 * L(k-1);
end
end