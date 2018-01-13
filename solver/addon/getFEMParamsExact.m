function [ matrix, eigenvalue, determinant ] = getFEMParamsExact(N)
switch N % for N<=4 exact values (rational) are exported from Wolfram Mathematica, help to ToMatlab.m package
    case 1
        A = [(1/2),(-1/2);(-1/2),(1/2)];
        C = [(2/3),(1/3);(1/3),(2/3)];
        eigenvalue.G = [0,6];
        determinant.C = (1/3);
        eigenvalue.Gov = [];
        determinant.Cov = 1;
    case 2
        A = [(7/6),(-4/3),(1/6);(-4/3),(8/3),(-4/3);(1/6),(-4/3),(7/6)];
        C = [(4/15),(2/15),(-1/15);(2/15),(16/15),(2/15);(-1/15),(2/15),(4/15) ...
            ];
        eigenvalue.G = [0,6,30];
        determinant.C = (8/135);
        eigenvalue.Gov = [5];
        determinant.Cov = (16/15);
    case 3
        A = [(37/20),(-189/80),(27/40),(-13/80);(-189/80),(27/5),(-297/80),( ...
            27/40);(27/40),(-297/80),(27/5),(-189/80);(-13/80),(27/40),( ...
            -189/80),(37/20)];
        C = [(16/105),(33/280),(-3/70),(19/840);(33/280),(27/35),(-27/280),( ...
            -3/70);(-3/70),(-27/280),(27/35),(33/280);(19/840),(-3/70),( ...
            33/280),(16/105)];
        eigenvalue.G = [0,45+(-1).*1605.^(1/2),30,45+1605.^(1/2)];
        determinant.C = (2187/224000);
        eigenvalue.Gov = [5,21];
        determinant.Cov = (6561/11200);
    case 4
        A = [(985/378),(-3424/945),(508/315),(-736/945),(347/1890);(-3424/945) ...
            ,(1664/189),(-2368/315),(2944/945),(-736/945);(508/315),( ...
            -2368/315),(248/21),(-2368/315),(508/315);(-736/945),(2944/945),( ...
            -2368/315),(1664/189),(-3424/945);(347/1890),(-736/945),(508/315), ...
            (-3424/945),(985/378)];
        C = [(292/2835),(296/2835),(-58/945),(8/405),(-29/2835);(296/2835),( ...
            256/405),(-128/945),(256/2835),(8/405);(-58/945),(-128/945),( ...
            208/315),(-128/945),(-58/945);(8/405),(256/2835),(-128/945),( ...
            256/405),(296/2835);(-29/2835),(8/405),(-58/945),(296/2835),( ...
            292/2835)];
        eigenvalue.G = [0,45+(-1).*1605.^(1/2),105+(-3).*805.^(1/2),45+1605.^(1/2),105+ ...
            3.*805.^(1/2)];
        determinant.C = (33554432/21097715625);
        eigenvalue.Gov = [28+(-2).*133.^(1/2),21,28+2.*133.^(1/2)];
        determinant.Cov = (67108864/281302875);
    otherwise
        [x, w ] = calcGaussianQuadrature(N+1);
        [f, df] = calcFiniteElementValue(x);
        C = NaN * ones(N+1); A = C;
        for i=1:(N+1) % SPEED UP: really, we need only 1/4 of calculated values
            for j=1:(N+1)
                C(i, j) = sum( f(:, i) .*  f(:, j) .* w);
                A(i, j) = sum(df(:, i) .* df(:, j) .* w);
            end
        end
end
matrix.A = A;
matrix.C = C;

%% Calc eigenvalues and determinants
if nargout >= 2 && N > 4
    if nargout == 3
        [ eigenvalue.G   determinant.C   ] = calc_eig(A, C);
        [ eigenvalue.Gov determinant.Cov ] = calc_eig(A(2:N, 2:N  ), C(2:N, 2:N  ));
    else % nargout == 2
        eigenvalue.G   = calc_eig(A, C);
        eigenvalue.Gov = calc_eig(A(2:N, 2:N  ), C(2:N, 2:N  ));
    end
end
end

function [G, C] = calc_eig(A, C)
[~, D] = eig(A, C/2, 'chol');
G = transpose(sort(diag(D)));
if nargout == 2
    C = det(C);
end
end