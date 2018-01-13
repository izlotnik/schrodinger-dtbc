function z = GaussianWave(pTask, Params, t, x, y)
switch pTask.sDimension
    case '1D'
        switch pTask.sExampleName
            case 'EX03'
                z = gw(x, t, Params) .* exp(1i * Params.exLambda .* t);
            case {'EX01r', 'EX20', 'EX20r', 'EX22r', 'EX24', 'EX25'} % normalized Gaussian beam
                z = gw(x, t, Params) ./ sqrt(sqrt(2 * pi * Params.exAlpha));
            otherwise
                z = gw(x, t, Params);
        end
    case '2D'
        x0 = Params.exX0; k_x = Params.exKoef ; alpha_x = Params.exAlpha ;
        y0 = Params.exY0; k_y = Params.exKoefY; alpha_y = Params.exAlphaY;
        [ X1, X2 ] = meshgrid(x, y);
        %
        pTask.sDimension = '1D'; p.exX0 = 0; 
        p.exKoef = k_x; p.exAlpha = alpha_x;
        z = bsxfun(@(x, t) GaussianWave(pTask, p, t, x), ( X1 - x0 + ( X2 - y0 ) ) / sqrt(2), t);
        %
        p.exKoef = k_y; p.exAlpha = alpha_y;
        z = bsxfun(@times, z, bsxfun(@(x, t) GaussianWave(pTask, p, t, x), ( X1 - x0 - ( X2 - y0 ) ) / sqrt(2), t));
end
end

function z = gw(x, t, Params)
x0    = Params.exX0;
k     = Params.exKoef;
alpha = Params.exAlpha;
%%
z = exp(1i .* k .* ( x - x0 - k .* t ) ...
    - ( x - x0 - 2 .* k .* t ) .^ 2 ./ ( 4 .* ( alpha + 1i .* t ) ) ...
    ) ./ sqrt(1 + 1i .* t ./ alpha);
end