function b = calcDST(a)
%% calcDST  
%
% function is based on MATLAB's dst function
[n, m] = size(a);
n = n + 1;
y = zeros(2*n, m);
y(  2:   n, :) =   a;
y(n+2: 2*n, :) = - a(end:-1:1,:); % - flipud(a);
z = fft(y);
b = z(2:n, :) / ( - 2 * 1i );