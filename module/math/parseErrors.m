function onorm = parseErrors(rError, params)
%%
error = containers.Map({'C', 'L2'}, { ...
    containers.Map({'abs', 'rel', 'sol'}, {rError.CAbs, rError.CRel, rError.CNormOfSolution}), ...
    containers.Map({'abs', 'rel', 'sol'}, {rError.L2Abs, rError.L2Rel, rError.L2NormOfSolution}) });
T = params('T');
M = params('M');
P = params('P');
L = size(P, 1);
onorm = containers.Map('KeyType', 'char', 'ValueType', 'any'); ;
for norm = error.keys()
    cnorm = error(norm{:});
    otype = containers.Map('KeyType', 'char', 'ValueType', 'any'); ;
    for type = cnorm.keys()
        ctype = cnorm(type{:});
        m = NaN(L, 1); t = m;
        for k = 1:L
            % One line MATLAB implementation of norm('L2')('abs'):
            % subsref(norm, struct('type', '()', 'subs', {'L2', {'abs'}}));
            [m(k), i] = parseError(ctype, fix(M*P(k, :)/T) + 1);
            t(k) = (i-1)/M*T;
        end
        otype(type{:}) = containers.Map({'error', 'time'}, {m, t});
    end
    onorm(norm{:}) = otype;
end

end

function [m, i] = parseError(error, j)
[m, i] = max(error(j(1):j(2))); i = j(1) + ( i - 1 );
end