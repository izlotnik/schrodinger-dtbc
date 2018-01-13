function value = GetValue(params, name, defaultValue)
if params.isKey(name)
    value = params(name);
else
    value = defaultValue; % return default value if unset
end;