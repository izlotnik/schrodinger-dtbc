function R = calcR(m, param1, param2)
if nargin <= 1
    error(message('calcR:TooFewInputs'));
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
                error(message('calcR:InvalidData'));
        end
    else if nargin == 3
            kappa = param1;
            mu    = param2;
        else
            error(message('calcR:TooMuchInputs'));
        end
    end
end
R = NaN * ones(m, 1);
R(1) = 1;
R(2) = - kappa * mu;
for k=2:(m-1)
    R(k+1) = ( ( 2 * k - 3 ) / k ) * kappa * mu * R(k) - ( ( k - 3 ) / k ) * kappa^2 * R(k-1);
end
end