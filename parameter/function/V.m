function z = V(pTask, x, y) % pTask.V_ -> pPotentialParams
switch pTask.sExampleName % sType : DBSQW / QB / PT / ST / ZERO
    case {'EX20','EX20r'}
        z = (                      x < pTask.V_x(1)) .* pTask.V_Q(1) ...
            + (x >= pTask.V_x(1) & x < pTask.V_x(2)) .* pTask.V_Q(2) ...
            + (x >= pTask.V_x(2) & x < pTask.V_x(3)) .* pTask.V_Q(3) ...
            + (x >= pTask.V_x(3) & x < pTask.V_x(4)) .* pTask.V_Q(4) ...
            + (x >= pTask.V_x(4) & x < pTask.V_x(5)) .* pTask.V_Q(5) ...
            + (x >= pTask.V_x(5)                   ) .* pTask.V_Q(6);
    case {'EX21', 'EX22', 'EX22a', 'EX22b', 'EX22m', 'EX22n', 'EX22k', 'EX22r'}
        z = (                      x < pTask.V_x(1)) .* pTask.V_Q(1) ...
            + (x >= pTask.V_x(1) & x < pTask.V_x(2)) .* pTask.V_Q(2) ...
            + (x >= pTask.V_x(2)                   ) .* pTask.V_Q(3);
    % case 'EX23' % delta-function
    case {'EX22s', 'EX24', 'EX25'}
        z = (                      x < pTask.V_x(1)) .* pTask.V_Q(1) ...
            + (x >= pTask.V_x(1)                   ) .* pTask.V_Q(2);
    case 'EX30'
        z = - pt(x, pTask);
    case 'EX31'
        z = st(x, pTask);
    case 'EX32'
        z = pw(x, pTask);
    otherwise
        z = repmat(pTask.V_Rinf, size(x));
end
if strcmp(pTask.sDimension, '2D') % add in 2D case
    z = z + 0; % instead of sTask.V_Q = pTask.V_Q + S(ll);
end
end

function z = pt(x, Params)
x0     = Params.V_x0;
alpha  = Params.V_alpha;
lambda = Params.V_lambda;
z = - alpha .^ 2 .* lambda .* (lambda + 1) .* sech( alpha .* (x - x0) ) .^ 2;
end

function z = st(x, Params)
x0     = Params.V_x0;
alpha  = Params.V_alpha;
z = ( Params.V_Rinf + Params.V_Linf ) / 2 + ( Params.V_Rinf - Params.V_Linf ) / 2 * tanh( alpha .* (x - x0) );
end

function z = pw(x, Params)
Q      = Params.V_Q; % 4000 % 1000
xL     = Params.V_x(1); % 1.9
xR     = Params.V_x(2); % 2.0
alpha  = Params.V_alpha; % 50
z = - Q * (tanh( alpha .* (x - xL) ) + 1) .* (tanh( alpha .* (xR - x) ) + 1);
end