function X = TRISYS(vLower, vDiag, vUpper, vRight)
%%
n = length(vRight);
% initialize
X = zeros(1, n); a = X; b = a;
% forward
g = vDiag(1);
a(1) = - vUpper(1) / g;
b(1) =   vRight(1) / g;
for k = 2:(n-1)
  g = vDiag(k) + vLower(k-1) * a(k-1);
  a(k) = - vUpper(k) / g;
  b(k) = ( vRight(k) - vLower(k-1) * b(k-1) ) / g;
end
% backward
X(n) = ( vRight(n) - vLower(n-1) * b(n-1) ) / (  vDiag(n) + vLower(n-1) * a(n-1) );
for k = (n-1):-1:1
  X(k) = a(k) * X(k+1) + b(k);
end