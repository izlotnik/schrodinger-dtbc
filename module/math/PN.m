function w = PN(D, n, h, x, p)
w = zeros(size(x));
for s=1:length(x);
    j = round(x(s)/h-1/2);
    if j < 0
        j = 0;
    elseif j >= p
        j = p - 1;
    end
    t = (x(s)-h*j)*n/h;
    for k=n:-1:0
        w(s) = w(s)*(t-k)/(k+1) + D(j+1, k+1);
    end
end
end