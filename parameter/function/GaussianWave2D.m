function z = GaussianWave2D( x, y, t, Params )
exX0     = Params.exX0;
exAlpha  = Params.exAlpha;
exKoef   = Params.exKoef;
%
exY0     = Params.exY0;
exAlphaY = Params.exAlphaY;
exKoefY  = Params.exKoefY;
%% EX01
% keyboard
% z = 1 / sqrt( 1 + 1i * t / exAlpha )...
%       * exp ( 1i * exKoef  .* ( x - exX0 - exKoef  * t )...
%       - ( x - exX0 - 2 * exKoef  * t ) .^ 2 / ( 4 * ( exAlpha + 1i * t ) ) )...
%       * conj(exp ( 1i * exKoefY .* ( y - exY0 - exKoefY * t )...
%       - ( y - exY0 - 2 * exKoefY * t ) .^ 2 / ( 4 * ( exAlphaY + 1i * t ) ) ))';
%% EX02
% z = 1 / sqrt( 1 + 1i * t / exAlpha ) ...
%       * exp ( 1i * exKoef .* ( x(:,1) - exX0 - exKoef * t )...
%       - ( x(:,1) - exX0 - 2 * exKoef * t ) .^ 2 / ( 4 * ( exAlpha + 1i * t ) ) )...
%       * conj( exp ( 1i * exKoefY .* ( x(:,2) - exY0 - exKoefY * t )...
%       - ( x(:,2) - exY0 - 2 * exKoefY * t ) .^ 2 / ( 4 * ( exAlphaY + 1i * t ) ) ) )';
% %% 
% % z = GaussianWave( x(:,1), t, p );
% %% 
% p.exX0     = 0;
% p.exAlpha  = Params.exAlphaY;
% p.exKoef   = Params.exKoefY;
% z = z * conj( GaussianWave( ( x(:,1) - exX0 + x(:,2) - exY0 ) / sqrt(2), t, p ) )';
% z = z * conj(GaussianWave( x(:,2), t, p ))';
%% EX02m
% % % OLD, TEST, 2010/02/02
% % % N = length(x); M = length(y);
% % % up = zeros(N,M); um = up;
% % % for n=1:N
% % %     for m=1:M
% % %         up(n,m) = ( x(n) - exX0 + ( y(m) - exY0 ) ) / sqrt(2);
% % %         um(n,m) = ( x(n) - exX0 - ( y(m) - exY0 ) ) / sqrt(2);
% % %     end
% % % end
% % % %%
% % % keyboard
% % % up = bsxfun( @plus , x , y' );
% % % um = bsxfun( @minus, x', y  );
%% EX02, EX22
[X1, X2] = meshgrid(x, y);
%%
p.exX0     = 0;
p.exAlpha  = exAlpha;
p.exKoef   = exKoef;
z = bsxfun(@(x,t) GaussianWave(x,t,p), ( X1 - exX0 + ( X2 - exY0 ) ) / sqrt(2), t);
%%
p.exX0     = 0;
p.exAlpha  = exAlphaY;
p.exKoef   = exKoefY;
z = bsxfun(@times,z,bsxfun(@(x,t) GaussianWave(x,t,p), ( X1 - exX0 - ( X2 - exY0 ) ) / sqrt(2), t) );
% %% EX03 (2D)
% z = bsxfun(@times,z,exp( 1i * Params.exLambda .* t ));
%% EX04
% [X1, X2] = meshgrid(x, y);
% %%
% p.exX0     = exX0;
% p.exAlpha  = exAlpha;
% p.exKoef   = exKoef;
% z = bsxfun(@(x,t) GaussianWave(x,t,p), X1, t);
% %%
% p.exX0     = exY0;
% p.exAlpha  = exAlphaY;
% p.exKoef   = exKoefY;
% z = bsxfun(@times,z,bsxfun(@(x,t) GaussianWave(x,t,p), X2, t) );