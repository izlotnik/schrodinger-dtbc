function params = prepareSolutionParams(params)
%
if params.isKey('N')
    params('matrix') = 'QR';
elseif params.isKey('theta')
    params('matrix') = 'TRISYS';
else
    params('matrix') = '*';
end

%
if params.isKey('r')
    params('method') = 'RICH';
else
    params('method') = 'CN';
end
end