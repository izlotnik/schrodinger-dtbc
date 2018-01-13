function [ K, c, R, L, s ] = calcConvolutionKernel(pTask, pCalc)
%%
hbar    = pTask.hbar;
B_inf   = pTask.B_inf;
rho_inf = pTask.rho_inf;
V_inf   = pTask.V_inf;
%%
h   = pCalc.h;
tau = pCalc.tau;
m   = pCalc.m;

%% Define family parameter
switch pCalc.name
    case {'FFDS', 'FFDS_M', 'FFDS_N', 'FFDS_Z'}
        theta = pCalc.theta;
    case {'FEM', 'FEM_C', 'FEM_F', 'FEM_Z'}
        N = pCalc.N;
end

%% Define a
switch pCalc.name
    case {'FFDS','FEM','OFF','FEM_C','FEM_F','FEM_Z'}
        a0 = (     h^2 * V_inf  ) / (     hbar^2 * B_inf );
        a1 = ( 2 * h^2 * rho_inf) / ( tau * hbar * B_inf );
        a  = a0 + 1i * a1;
    case {'FFDS_M','FFDS_N'}
        a0 = (       V_inf ) / (     hbar^2 * B_inf );
        a1 = ( 2 * rho_inf ) / ( tau * hbar * B_inf );
        a  = h^2*(a0 + 1i * a1);
    case {'SD','FFDS_Z'}
        a0 = (       V_inf ) / (     hbar^2 * B_inf );
        a1 = ( 2 * rho_inf ) / ( tau * hbar * B_inf );
        a  = a0 + a1 * 1i;
end

%% Calc kernel
switch pCalc.name
    case 'FFDS_Z'
        D = [ a 2+(1-4*theta)*h^2*a ];
        K = calcR(m, arg(D./conj(D)));
        if nargout >= 2
            A = D(1)*D(2);
            B = D(1)/D(2);
            if imag(B) > 0
                s = -1;
            else
                s =  1;
            end
            c = s * sqrt(abs(A))*exp(-1i*arg(A)/2)/2;
            if nargout >= 3
                R = K;
                if nargout >= 4
                    L = [];
                end
            end
        end
    case 'FFDS_M'
        D = [ a 2+(1-4*theta)*a ];
        K = calcR(m, arg(D./conj(D)));
        if nargout >= 2
            A = conj(D(1)*D(2));
            B = conj(D(1)/D(2));
            if imag(B) < 0
                s = -1;
            else
                s =  1;
            end
            c = - s * sqrt(abs(A))*exp(1i*arg(A)/2)/(2*h);
            if nargout >= 3
                R = K;
                if nargout >= 4
                    L = [];
                end
            end
        end
    case 'FFDS_N'
        D = [ a 2+(1-4*theta)*a ];
        K = calcR(m, arg(D./conj(D)));
        if nargout >= 2
            A = conj(prod(D));
            B = sqrt(abs(A))*exp(1i*arg(A)/2);
            g = conj([ (D(1)+D(2))/2 (D(1)-D(2))/2 ]); % g = conj([ 1+(1-2*theta)*a 2*theta*a-1 ]);
            if abs(-g(1)-B) < abs(g(2))
                s = -1;
            else
                if abs(-g(1)+B) < abs(g(2))
                    s = 1;
                else
                    s = NaN; warning(message('calcKernel:UnexpectedValue'));
                end
            end
            c = - s * B/(2*h);
            if nargout >= 3
                R = K;
                if nargout >= 4
                    L = [];
                end
            end
        end
    case 'FFDS'
        b     = 1 - 2 * theta * a;
        alpha = a * ( a + 2 * b );
        beta  = 2 * a0 + ( 1 - 4 * theta ) * abs(a)^2;
        arga  = arg(alpha);
        mu    = beta / abs(alpha);
        if abs(mu) >= 1, warning(message('calcKernel:UnexpectedValue')); end
        kappa = - exp(1i*arga);
        K = calcR(m, kappa, mu);
        if nargout >= 2
            argb  = arg(b);
            if 2*argb - arga < 0
                k0 = -1;
            else
                if 2*argb - arga < 2*pi
                    k0 = 0;
                else
                    if 2*argb - arga < 4*pi
                        k0 = 1;
                    else
                        k0 = NaN; warning(message('calcKernel:UnexpectedValue'));
                    end
                end
            end
            s = (-1)^k0;
            c = s * sqrt(abs(alpha))*exp(-1i*arga/2)/(2*h);
            if nargout >= 3
                R = K;
                if nargout >= 4
                    L = [];
                end
            end
        end
    case 'OFF'
        alpha = a*(2+a);
        beta  = 2*a0+abs(a)^2;
        arga  = arg(alpha);
        mu    = beta / abs(alpha);
        if abs(mu) >= 1, warning(message('calcKernel:UnexpectedValue')); end
        kappa = - exp(1i*arga);
        K = calcR(m, kappa, mu);
        if nargout >= 2
            s = -1;
            c = s * sqrt(abs(alpha))*exp(-1i*arga/2)/(2*h);
            if nargout >= 3
                R = K;
                if nargout >= 4
                    L = [];
                end
            end
        end
    case 'SD'
        arga = arg(a);
        mu   = a0 / abs(a);
        if abs(mu) >= 1, warning(message('calcKernel:UnexpectedValue')); end
        kappa = - exp(1i*arga);
        K = calcR(m, kappa, mu);
        if nargout >= 2
            s = -1;
            c = s * sqrt(abs(a)/2)*exp(-1i*arga/2);
            if nargout >= 3
                R = K;
                if nargout >= 4
                    L = [];
                end
            end
        end
    case 'FEM'
        [ M, E ] = getFEMParams(N);
        n1 = fix((N+1)/2);
        R = NaN * ones(m, n1); L = R;
        s = n1;
        R(:, s) = calcR(m, arg((a+E.G((2*s-1):2*s)) ./ (conj(a)+E.G((2*s-1):2*s))));
        K = R(:, s);
        if mod(N, 2) == 0
            E.Gov(2*s) = E.G(2*s+1);
            L(:, s) = calcL(m, arg((a+E.Gov((2*s-1):2*s)) ./ (conj(a)+E.Gov((2*s-1):2*s))));
            Q = zeros(m, 1); Q(1) = 1; Q(2) = exp(1i*arg((a+E.G(2*s+1))/(conj(a)+E.G(2*s+1))));
            Q = conv2(L(:, s), Q);              % Q = conv_fft2(L(:, s), Q);
            K = conv2(K, Q(1:m)); K = K(1:m);   % K = conv_fft2(K, Q(1:m));
        end
        for s = 1:(n1-1)
            R(:, s) = calcR(m, arg((a+E.G((2*s-1):2*s)) ./ (conj(a)+E.G((2*s-1):2*s))));
            K = conv2(K, R(:, s)); K = K(1:m);  % K = conv_fft2(K, R(:, s));
            L(:, s) = calcL(m, arg((a+E.Gov((2*s-1):2*s)) ./ (conj(a)+E.Gov((2*s-1):2*s))));
            K = conv2(K, L(:, s)); K = K(1:m);  % K = conv_fft2(K, L(:, s));
        end
        if nargout >= 2
            if N > 1
                Gov = det(M.A(2:N, 2:N) + conj(a)/2 * M.C(2:N, 2:N));
                g = [ det(M.A(1:N, 1:N) + conj(a)/2 * M.C(1:N, 1:N)) ...
                    (-1)^(N-1)*det(M.A(1:N, 2:(N+1)) + conj(a)/2 * M.C(1:N, 2:(N+1))) ];
            else
                Gov = 1;
                g = [ M.A(1, 1) + conj(a)/2 * M.C(1, 1) ...
                    M.A(1, 2) + conj(a)/2 * M.C(1, 2) ];
            end
            g = g / Gov;
            A = det(M.A(:, :) + conj(a)/2 * M.C(:, :))/Gov;
            B = sqrt(abs(A))*exp(1i*arg(A)/2);
            if abs(-g(1)-B) < abs(g(2))
                s = -1; %ALWAYS%
            else
                if abs(-g(1)+B) < abs(g(2))
                    s = 1;
                else
                    s = NaN; warning(message('calcKernel:UnexpectedValue'));
                end
            end
            c = - s*B/h;
        end
    case 'FEM_C'
        [ M, E ] = getFEMParams(N);
        n1 = fix((N+1)/2);
        R = NaN * ones(m, n1); L = R;
        s = n1;
        R(:, s) = calcR(m, arg((a+E.G((2*s-1):2*s)) ./ (conj(a)+E.G((2*s-1):2*s))));
        K = R(:, s);
        if mod(N, 2) == 0 % N - четное
            E.Gov(2*s) = E.G(2*s+1);
            % L(:, s) = calcL(m, arg([ (a+E.Gov(2*s-1)) / (conj(a)+E.Gov(2*s-1)) (a+E.Gov(2*s)) / (conj(a)+E.Gov(2*s)) ]));
            L(:, s) = calcL(m, arg((a+E.Gov((2*s-1):2*s)) ./ (conj(a)+E.Gov((2*s-1):2*s))));
            Q = zeros(m, 1); Q(1) = 1; Q(2) = exp(1i*arg((a+E.G(2*s+1))/(conj(a)+E.G(2*s+1))));
            Q = conv2(L(:, s), Q);
            K = conv2(K, Q(1:m)); K = K(1:m);
            % L2 = NaN * ones(m, 1); L2(1) = L(1, s);
            % L2(2:m) = L(2:m, s)+(a+E.Gov(2*s))/(conj(a)+E.Gov(2*s))*L(1:(m-1),s);
            % K  = conv_fft2(K, L2);
        end
        for s = 1:(n1-1)
            R(:, s) = calcR(m, arg((a+E.G((2*s-1):2*s)) ./ (conj(a)+E.G((2*s-1):2*s))));
            % R(:, s) = calcR(m, arg([ (a+E.G(2*s-1)) / (conj(a)+E.G(2*s-1)) (a+E.G(2*s)) / (conj(a)+E.G(2*s)) ]));
            K = conv2(K, R(:, s)); K = K(1:m);
            L(:, s) = calcL(m, arg((a+E.Gov((2*s-1):2*s)) ./ (conj(a)+E.Gov((2*s-1):2*s))));
            % L(:, s) = calcL(m, arg([ (a+E.Gov(2*s-1)) / (conj(a)+E.Gov(2*s-1)) (a+E.Gov(2*s)) / (conj(a)+E.Gov(2*s)) ]));
            K = conv2(K, L(:, s)); K = K(1:m);
        end
        if nargout >= 2
            if N > 1
                Gov = det(M.A(2:N, 2:N) + conj(a)/2 * M.C(2:N, 2:N));
                g = [ det(M.A(1:N, 1:N) + conj(a)/2 * M.C(1:N, 1:N)) ...
                    (-1)^(N-1)*det(M.A(1:N, 2:(N+1)) + conj(a)/2 * M.C(1:N, 2:(N+1))) ];
            else
                Gov = 1;
                g = [ M.A(1, 1) + conj(a)/2 * M.C(1, 1) ...
                    M.A(1, 2) + conj(a)/2 * M.C(1, 2) ];
            end
            g = g / Gov;
            A = det(M.A(:, :) + conj(a)/2 * M.C(:, :))/Gov;
            B = sqrt(abs(A))*exp(1i*arg(A)/2);
            if abs(-g(1)-B) < abs(g(2))
                s = -1; %ALWAYS%
            else
                if abs(-g(1)+B) < abs(g(2))
                    s = 1;
                else
                    s = NaN; warning(message('calcKernel:UnexpectedValue'));
                end
            end
            c = - s*B/h;
        end
    case 'FEM_F'
        [ M, E ] = getFEMParams(N);
        n1 = fix((N+1)/2);
        R = NaN * ones(m, n1); L = R;
        s = n1;
        R(:, s) = calcR(m, arg((a+E.G((2*s-1):2*s)) ./ (conj(a)+E.G((2*s-1):2*s))));
        K = R(:, s);
        if mod(N, 2) == 0 % N - четное
            E.Gov(2*s) = E.G(2*s+1);
            L(:, s) = calcL(m, arg((a+E.Gov((2*s-1):2*s)) ./ (conj(a)+E.Gov((2*s-1):2*s))));
            Q = zeros(m, 1); Q(1) = 1; Q(2) = exp(1i*arg((a+E.G(2*s+1))/(conj(a)+E.G(2*s+1))));
            Q = conv_fft2(L(:, s), Q);
            K = conv_fft2(K, Q(1:m)); K = K(1:m);
        end
        for s = 1:(n1-1)
            R(:, s) = calcR(m, arg((a+E.G((2*s-1):2*s)) ./ (conj(a)+E.G((2*s-1):2*s))));
            K = conv_fft2(K, R(:, s)); K = K(1:m);
            L(:, s) = calcL(m, arg((a+E.Gov((2*s-1):2*s)) ./ (conj(a)+E.Gov((2*s-1):2*s))));
            K = conv_fft2(K, L(:, s)); K = K(1:m);
        end
        if nargout >= 2
            if N > 1
                Gov = det(M.A(2:N, 2:N) + conj(a)/2 * M.C(2:N, 2:N));
                g = [ det(M.A(1:N, 1:N) + conj(a)/2 * M.C(1:N, 1:N)) ...
                    (-1)^(N-1)*det(M.A(1:N, 2:(N+1)) + conj(a)/2 * M.C(1:N, 2:(N+1))) ];
            else
                Gov = 1;
                g = [ M.A(1, 1) + conj(a)/2 * M.C(1, 1) ...
                    M.A(1, 2) + conj(a)/2 * M.C(1, 2) ];
            end
            g = g / Gov;
            A = det(M.A(:, :) + conj(a)/2 * M.C(:, :))/Gov;
            B = sqrt(abs(A))*exp(1i*arg(A)/2);
            if abs(-g(1)-B) < abs(g(2))
                s = -1; %ALWAYS%
            else
                if abs(-g(1)+B) < abs(g(2))
                    s = 1;
                else
                    s = NaN; warning(message('calcKernel:UnexpectedValue'));
                end
            end
            c = - s*B/h;
        end
end
end