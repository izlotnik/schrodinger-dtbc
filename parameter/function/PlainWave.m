function z = PlainWave(pTask, Params, t, x, y)
switch pTask.sDimension
    case '1D'
        exX0    = Params.exX0;  exKoef  = Params.exKoef;
        exBeta  = exKoef;       exAlpha = Params.exAlpha;
        %%
        z = exp(1i .* exKoef .* ( x - exX0 - exBeta .* t )) * exAlpha ;
    case '2D'
        x0 = Params.exX0; k_x = Params.exKoef ; alpha_x = Params.exAlpha ;
        y0 = Params.exY0; k_y = Params.exKoefY; alpha_y = Params.exAlphaY;
        [ X1, X2 ] = meshgrid(x, y);
        %
        pTask.sDimension = '1D'; p.exX0 = 0;
        p.exKoef = k_x; p.exAlpha = alpha_x;
        z = bsxfun(@(x, t) PlainWave(pTask, p, t, x), ( X1 - x0 + ( X2 - y0 ) ) / sqrt(2), t);
        %
        p.exKoef = k_y; p.exAlpha = alpha_y;
        z = bsxfun(@times, z, bsxfun(@(x, t) PlainWave(pTask, p, t, x), ( X1 - x0 - ( X2 - y0 ) ) / sqrt(2), t));
end
end