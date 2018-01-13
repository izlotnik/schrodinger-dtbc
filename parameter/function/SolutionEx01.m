function z = SolutionEx01( x, y, t )
% exX0     = Params.exX0;
% exAlpha  = Params.exAlpha;
% exKoef   = Params.exKoef;
% %
% exY0     = Params.exY0;
% exAlphaY = Params.exAlphaY;
% exKoefY  = Params.exKoefY;
%%
z = 1i * exp(1i*t) ./ ( cosh(x).* cosh(y) );