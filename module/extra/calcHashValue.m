function hash = calcHashValue(pCalc, pTask)
hash = str2hash([getStr('', pCalc) '#' getStr('', pTask)]);
end

function str = getStr(str, val)
sSeparator = '@';
if islogical(val) || ( isscalar(val) && isinteger(val) )
    str = [ str num2str(val) ];
elseif ischar(val)
    str = [ str val ];
elseif isscalar(val) && isfloat(val)
    str = [ str num2str(val, '%20.6g') ];
elseif isstruct(val)
    sFields = sort(fieldnames(val));
    sFields = sFields(~strcmp(sFields, 'sBaseRunID'));
    sFields = sFields(~strcmp(sFields, 'sPath'));
    for k=1:length(sFields)
        str = [ getStr(str, getfield(val, sFields{k})) sSeparator ];
    end
elseif isvector(val)
    for k=1:length(val(:))
        str = [ str num2str(val(k), '%20.6g') sSeparator ];
    end
end
end

function str = str2hash(str)
% str = MessageDigest(str);
% str = hash(str, 'MD5'); 
str = DataHash(str, ...
    struct('Format', 'hex', 'Method', 'MD5', 'Input', 'ascii'));
end