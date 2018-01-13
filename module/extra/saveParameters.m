function saveParameters(sFileName, Params)
sFields = sort(fieldnames(Params));
fileID = fopen(sFileName, 'a');
for k=1:length(sFields)
    fieldValue = getfield(Params, sFields{k});
    if islogical(fieldValue) || ( isscalar(fieldValue) && isinteger(fieldValue) )
        fprintf(fileID, '%s=' , sFields{k});
        fprintf(fileID, '%d\n', fieldValue);
    elseif ischar(fieldValue)
        fprintf(fileID, '%s=' , sFields{k});
        fprintf(fileID, '%s\n', fieldValue);
    elseif isscalar(fieldValue) && isfloat(fieldValue)
        fprintf(fileID, '%s=' , sFields{k});
        fprintf(fileID, '%s\n', num2str(fieldValue,'%20.20g'));
    end
end
fprintf(fileID, '%s\n', '');
fclose(fileID);