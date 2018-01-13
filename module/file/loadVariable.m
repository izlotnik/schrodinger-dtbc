function out = loadVariable(sFile, bFlag, Params)
% sLogType      = 'info'; sLogMessageID = 'LoadData';
if bFlag
    % logMessage(sLogType, sLogMessageID, 'Loading input data in %s'         , sFile, Params);
    out = load(sFile); 
    if isfield(out, 'pVariable')
        out = out.pVariable;
    end
    % logMessage(sLogType, sLogMessageID, 'Input data is loaded from file %s', sFile, Params);
else
    out = [];
end