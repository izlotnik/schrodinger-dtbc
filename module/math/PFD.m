function D = PFD(y, p, n)
D = NaN * ones(p, n+1);
for r=0:(p-1)
    D(r+1, :) = FD(y((r*n+1):((r+1)*n+1)), n);
end
end

function d = FD(y, n)
d = NaN * ones(n+1, n+1);
d(:, 1) = y(:);
for k=1:n
    for i=0:(n-k)
        d(i+1, k+1) = d(i+2, k) - d(i+1, k);
    end
end
d = d(1, :);
end