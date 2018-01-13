function logMessage(sType, sMessageID, sMessage, value, sMFunctionName)
sProject = 'dtbc';
%% message ID must have the form: '<sProject>:<sMFunctionName>:<sMessageID>'
if isempty(sMFunctionName)
    sMessageIDOut  = sMessageID;
    sMessagePrefix = [];
else
    sMessageIDOut  = [ sMFunctionName ':' sMessageID ];
    sMessagePrefix = [ '[' sMFunctionName '.m]: ' ];
end
if isempty(findstr([ sProject ':' ], sMessageIDOut))
    sMessageIDOut  = [ sProject ':' sMessageIDOut ];
end
%% Form message in the form: '[<sMFunctionName>.m]: <sMessage>'
sMessageOut = sprintf([ sMessagePrefix sMessage ], value);
switch lower(sType)
    case {'info', 'ok', 'i'}
        disp([ 'Informational: ' sMessageOut ]);
    case {'warning', 'warn', 'w'}
        s = warning('query', 'backtrace');
        if strcmp(s.state,'on')
            warning('off', 'backtrace');
        end
        warning(sMessageIDOut, sMessageOut); warning(s);
    case {'error', 'err', 'e'}
        error(sMessageIDOut, sMessageOut);
    otherwise
        disp([ 'Informational: ' sMessageOut ]);
        logMessage('warning', sMessageID, 'Unknown message type.\n Please use one of the predefined type: ''info'',''warning'' or ''error''.', [], mfilename)
end