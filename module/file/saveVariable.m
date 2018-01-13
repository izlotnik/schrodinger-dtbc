function saveVariable(pVariable, sFile, bFlag, Params)
% sLogType      = 'info'; sLogMessageID = 'SaveData';
if bFlag
    % logMessage(sLogType, sLogMessageID, 'Saving output data in %s'         , sFile, Params);
    save(sFile, 'pVariable', '-v7.3');
    % logMessage(sLogType, sLogMessageID, 'Output data is saved into file %s', sFile, Params);
end