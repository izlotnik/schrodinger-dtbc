function sFile = saveImage(hFigure, pParams, sMnemo, Params)
if pParams.bSaveFigure
    sLogType      = 'info';
    sLogMessageID = 'SaveFigure';
    set(hFigure, 'PaperPositionMode', 'auto');
    % orient landscape;
    if isempty(sMnemo)
        sFileSubSuffix = sMnemo;
    else
        sFileSubSuffix = [ '_' sMnemo ]; 
    end
    switch pParams.sPrintType
        case '-depsc2'
            sFileSubSuffix = [ sFileSubSuffix '_' 'C'  ];
        case '-deps2'
            sFileSubSuffix = [ sFileSubSuffix '_' 'BW' ];
    end
    sFile = [ pParams.sPath '/' pParams.sFilePrefix '_' pParams.sFileName '_' pParams.sFileSuffix sFileSubSuffix ];
    % logMessage(sLogType, sLogMessageID, 'Saving figure into file\n%s', [ sFile '.' pParams.sFileExtension ], Params);
    fprintf('Save: %s\n', [ sFile '.' pParams.sFileExtension ]);
    i = 1; % print safely figure in maxumum 10 attempts
    while ~printSafe(hFigure, pParams.sPrintType, sFile, pParams.sFileExtension) && i < 10 
        i = i + 1;
    end
    % logMessage(sLogType, sLogMessageID, 'Figure is saved', '', Params);
    %
    % X = imread([ sFile '.' 'png' ]);
    % imwrite(X, [ sFile '_LOSELESS.' 'jpg' ], 'Mode', 'lossless');
    % imwrite(X, [ sFile '_QUALITY_100.' 'jpg' ], 'Quality', 100);
end
end

function status = printSafe(hFigure, sPrintType, sFile, sFileExtension)
status = 1;
try   
    % saveas(hFigure, [ sFile '.jpg' ], 'jpg')
    % saveas(hFigure, [ sFile '.png' ], 'png')
    % saveas(hFigure, [ sFile '.pdf' ], 'pdf')
    print(hFigure, sPrintType, [ sFile '.' sFileExtension ]);
    % print(hFigure, '-dpng', [ sFile '.png' ]);
catch exception % disp(exception.identifier)
    if exist([ sFile '.' 'epstmp' ], 'file')
        delete([ sFile '.' 'epstmp' ])
        status = 0;
    end
end
end