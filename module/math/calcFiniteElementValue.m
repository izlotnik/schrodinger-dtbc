function [f, df, y] = calcFiniteElementValue(x, m)
%calcFiniteElementValue
%   Return function (f) and its derivative's values (df) at x.

%% Check input values
if nargin < 2, m = length(x); end

%% Calculate f and its derivative at x
y = linspace(-1, 1, m); % -1 + 2*(i-1)/m, 1<=i<=m+1
f    = ones (m, m);
df   = zeros(m, m);
znam = ones (m, 1);
for k=1:m
    ss = zeros(m, 1);
    for s=1:m
        g = ones(m, 1);
        if s ~= k
            f(:, k) = f(:, k) .* (    x - y(s) );
            znam(k) = znam(k)  * ( y(k) - y(s) );
            for j=1:m
                if ( j ~= k ) && ( j ~= s )
                    g = g .* ( x - y(j) );
                end
            end
            ss = ss + g;
        end
    end
    f(:, k) = f(:, k) ./ znam(k);
    df(:,k) = ss      ./ znam(k);
end