function d = getDataPath()
[~, hostname] = system('hostname');
switch hostname(1:(end-1))
    case 'asus-u36sg'
        d = 'D:\Matlab\R2011b';
    case 'PC-MM'
        d = 'E:\Matlab\R2011b';
    otherwise
        d = 'C:\Matlab\R2011b';
end
end