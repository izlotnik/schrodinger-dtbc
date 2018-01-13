function printMaxErrors(rError, rModulError, Params)
switch Params.name
    case 'FEM'
        param = Params.N; type = 'n=%d';
    case 'FFDS'
        param = Params.theta; type = 'theta=%.4g';
end
logMessage('info', 'MaxErrorOutput'...
    , sprintf('\nFor J=%d (h=%.4g), M=%d (tau=%.4g) and %s maximum error values are:'...
    , Params.n-1, Params.h, Params.m-1, Params.tau, type), param, mfilename);
if ~isempty(rError)
    fprintf('max E_{L2}    =%d\nmax E_{C}     =%d\nmax E_{L2,Rel}=%d\nmax E_{C, Rel}=%d\n', ...
        max(rError.L2Abs), max(rError.CAbs), max(rError.L2Rel), max(rError.CRel));
end
if ~isempty(rModulError)
    fprintf('max E_{L2}    =%d\nmax E_{C}     =%d\nmax E_{L2,Rel}=%d\nmax E_{C, Rel}=%d\n', ...
        max(rModulError.L2Abs), max(rModulError.CAbs), max(rModulError.L2Rel), max(rModulError.CRel));
end