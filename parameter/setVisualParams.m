function [ pVisual, pCalc ] = setVisualParams( pTask, pCalc, pSolution, pError, pAddon )
%% Plot solution
pVisual.bPlotTimeBehaviour = pSolution(1);
bPlotPotential             = pSolution(2);
pVisual.nPlotCount         = pSolution(3);
pVisual.bPlotBehaviour3D   = pSolution(4);
pVisual.bLogSolution       = pSolution(5);
%% Plot error
pVisual.bPlotErrorTimeBehaviour        = pError(1);
nFrequency                             = pError(2); 
if ~isnan(nFrequency) && nFrequency == int32(nFrequency) && nFrequency > 0 && nFrequency < pCalc.m
    pCalc.nFreq = 1:nFrequency:pCalc.m;
else
    pCalc.nFreq = 1:pCalc.m;
end
pVisual.bPlotMaxError                  = pError(3);
pVisual.bPlotErrorModulusTimeBehaviour = pError(4);
%%
pVisual.bPlotDTBCKernel = pAddon(1);
bSaveFigure             = pAddon(2);
%% Movie parameters
pVisual.bPlayMovie = false;
pVisual.bSaveMovie = false;
pVisual.pPlot2D.pMovie.bSave = pVisual.bSaveMovie;
pVisual.pPlot3D.pMovie.bSave = false;
pVisual.pPlotError.bSave = false;
%%
sSaveImagePath = [ cd '\data\plots' ];    % pCalc.sWorkDirectory; %
sBaseStyle     = 'phd-2012'; %'report-2010'; % 'master-2010';   % 'vestnik-mpei-2010'; % 'krm-2009'; %
if pVisual.bPlotErrorTimeBehaviour
    pCalc.bCalcError = true;
%     logMessage( 'info', 'ChangeBooleanParameter', [ 'We need plotting errors so we must calculate them.\n'...
%         ' Variable bCalcError is set to true.' ], [], mfilename );
end
if pVisual.bPlotErrorModulusTimeBehaviour
    pCalc.bCalcErrorModulus = true;
%     logMessage( 'info', 'ChangeBooleanParameter', [ 'We need plotting absolute errors so we must calculate them.\n'...
%         ' Variable bCalcErrorModulus is set to true.'], [], mfilename );
end
if pVisual.bPlotBehaviour3D && strcmp(pTask.sDimension, '1D')
    pCalc.bStoreAllTimeLevels = true;
%     logMessage('info', 'ChangeValue', 'We need plotting 3D solution, so we need to store all time level solution.', [], mfilename);
end
switch pTask.sDimension
    case '1D'
        switch pCalc.name
            case 'FFDS'
                thetaSymbolPrint = sym2str(sym(pCalc.theta));
                pVisual.thetaSymbolPrint = thetaSymbolPrint;
        end
    case '2D'
        switch pCalc.name
            case 'FEM' % TODO: change it
                etaSymbolPrint   = sym2str(sym(pCalc.eta));
                pVisual.etaSymbolPrint   = etaSymbolPrint;
            case 'FFDS'
                thetaSymbolPrint = sym2str(sym(pCalc.theta));
                etaSymbolPrint   = sym2str(sym(pCalc.eta));
                pVisual.thetaSymbolPrint = thetaSymbolPrint;
                pVisual.etaSymbolPrint   = etaSymbolPrint;
        end
end
if pVisual.bPlotDTBCKernel
    pVisual.pPlotKernel.nFontSize       = 24;
    pVisual.pPlotKernel.sFontWeight     = 'bold';
    pVisual.pPlotKernel.sLineType       = '-b';
    pVisual.pPlotKernel.nLineWidth      = 3;
    pVisual.pPlotKernel.sUnits          = 'normalized';
    pVisual.pPlotKernel.sNextPlot       = 'add';
    pVisual.pPlotKernel.sNumberTitle    = 'off';
    pVisual.pPlotKernel.sMenuBar        = 'none';
    pVisual.pPlotKernel.sName           = '2D';
    pVisual.pPlotKernel.nOuterPosition  = [0 0 1 1];
    pVisual.pPlotKernel.nPosition       = [0 0.04 1 0.935];
    pVisual.pPlotKernel.sLegendLocation = 'SouthWest';
    pVisual.pPlotKernel.sStyle          = sBaseStyle;
    pVisual.pPlotKernel.bSaveFigure     = bSaveFigure;
    pVisual.pPlotKernel.sPath           = sSaveImagePath;
    pVisual.pPlotKernel.sPrintType      = '-depsc2';
    pVisual.pPlotKernel.sFilePrefix     = [ datestr(now, 'yyyymmdd_HHMMSS') '_' pTask.sDimension '_' pTask.sExampleName ];
    pVisual.pPlotKernel.sFileName       = 'DTBC_Kernel';
    pVisual.pPlotKernel.sLabelX         = '$m$';
    switch pCalc.name
        case 'FEM'
            pVisual.pPlotKernel.sFileSuffix = [ pCalc.name '_N=' num2str(pCalc.N) ];
            pVisual.pPlotKernel.sLegendTitle= {[ '$\left|K^{(' num2str(pCalc.N) '),m}_{\rm{ref}}\right|$' ]};
            pVisual.pPlotKernel.sLabelY     = '$\left| K^{(n),m}_{\rm ref}\right|$';
        case 'FFDS'
            pVisual.pPlotKernel.sFileSuffix = [ 'PARAM=' strrep(thetaSymbolPrint,'/','-') ];
            pVisual.pPlotKernel.sLegendTitle= {['$\left|R^{m}_{\theta=' thetaSymbolPrint '}\right|$']};
            pVisual.pPlotKernel.sLabelY     = '$\left| R^{m}_{\theta}\right|$';
        otherwise
            logMessage('error', 'WrongCalcName', 'Unknown calculation name.\n Current value is %s', pCalc.name, mfilename);
    end
    if strcmp(pCalc.type,'SSP')
        pVisual.pPlotKernel.sFileSuffix = [ pCalc.type '_' pVisual.pPlotKernel.sFileSuffix ];
    end
    pVisual.pPlotKernel.sFileExtension  = 'eps';
end
if pVisual.bPlotMaxError
    pVisual.pPlotMaxError.nFontSize       = 24;
    pVisual.pPlotMaxError.sFontWeight     = 'bold';
    pVisual.pPlotMaxError.sLineType       = '-b';
    pVisual.pPlotMaxError.nLineWidth      = 3;
    pVisual.pPlotMaxError.sUnits          = 'normalized';
    pVisual.pPlotMaxError.sNextPlot       = 'add';
    pVisual.pPlotMaxError.sName           = '2D';
    pVisual.pPlotMaxError.nOuterPosition  = [0 0 1 1];
    pVisual.pPlotMaxError.sStyle          = sBaseStyle;
    pVisual.pPlotMaxError.bSaveFigure     = bSaveFigure;
    pVisual.pPlotMaxError.sPath           = sSaveImagePath;
    pVisual.pPlotMaxError.sPrintType      = '-depsc2';
    pVisual.pPlotMaxError.sFilePrefix     = [ datestr(now, 'yyyymmdd_HHMMSS') '_' pTask.sDimension '_' pTask.sExampleName ];
    pVisual.pPlotMaxError.sFileName       = 'MaxError';
    pVisual.pPlotMaxError.sLabelX         = '$J$';
    pVisual.pPlotMaxError.sLabelY         = {[ '$\lg\max\limits_{1\leq m\leq ' num2str(pCalc.m-1) '} E^{(n),m}_{\rm rel}(J)$ ']};
    switch pCalc.name
        case 'FEM'
            pVisual.pPlotMaxError.sLegendTitle= [ '$m=' num2str(pCalc.m-1) ', n=' num2str(pCalc.N) '$' ];
        case 'FFDS'
            pVisual.pPlotMaxError.sLegendTitle= [' $m=' num2str(pCalc.m-1) ', \theta=' thetaSymbolPrint '$' ];
        otherwise
            logMessage('error', 'WrongCalcName', 'Unknown calculation name.\n Current value is %s', pCalc.name, mfilename);
    end
    pVisual.pPlotMaxError.sFileExtension  = 'eps';
end
if pVisual.bPlotTimeBehaviour
    switch pTask.sDimension
        case '1D'
            pVisual.bPlotTimeBehaviour2D    = true;
            pVisual.pPlot2D.nPlotCount      = pVisual.nPlotCount;
            pVisual.pPlot2D.bPlayPause      = false;
            pVisual.pPlot2D.nFontSize       = 24;
            pVisual.pPlot2D.sFontWeight     = 'bold';
            pVisual.pPlot2D.sLineType       = {'--r', '-b'};
            pVisual.pPlot2D.nLineWidth      = 3;
            if strcmp(pTask.sExampleName, 'EX01r')
                pVisual.pPlot2D.nPlotX      = -pCalc.X_0/2:pCalc.h_mod:pCalc.X_0/2;
                pVisual.pPlot2D.nLimitX     = [-pCalc.X_0/2 pCalc.X_0/2];
                pVisual.pPlot2D.nLimitY     = [-2.15 2.15];
            elseif strcmp(pTask.sExampleName, 'EX20r')
                pVisual.pPlot2D.nPlotX      = -pCalc.X_0/2:pCalc.h_mod:pCalc.X_0/2;
                pVisual.pPlot2D.nPlotM      = pCalc.x_mod;
                pVisual.pPlot2D.nLimitX     = [-pCalc.X_0/2 pCalc.X_0/2];
                pVisual.pPlot2D.nLimitY     = [-1 1];
            elseif strcmp(pTask.sExampleName, 'EX22m')
                pVisual.pPlot2D.nPlotX      = -pCalc.X_0/2:pCalc.h_mod:pCalc.X_0/2;
                pVisual.pPlot2D.nPlotM      = pCalc.x_mod;
                pVisual.pPlot2D.nLimitX     = [-pCalc.X_0/2 pCalc.X_0/2];
                pVisual.pPlot2D.nLimitY     = [-1.15 1.15];
            elseif strcmp(pTask.sExampleName, 'EX22r')
                pVisual.pPlot2D.nPlotX      = -pCalc.X_0/2:pCalc.h_mod:pCalc.X_0/2;
                pVisual.pPlot2D.nPlotM      = pCalc.x_mod;
                pVisual.pPlot2D.nLimitX     = [-pCalc.X_0/2 pCalc.X_0/2];
                pVisual.pPlot2D.nLimitY     = [-2.15 2.15];
            else
                pVisual.pPlot2D.nPlotX      = pCalc.x_mod; % pCalc.xInner; %
                pVisual.pPlot2D.nLimitX     = [0 pCalc.X_0];
                pVisual.pPlot2D.nLimitY     = [-1 1];
            end
            pVisual.pPlot2D.sUnits          = 'normalized';
            pVisual.pPlot2D.sNextPlot       = 'replacechildren';
            pVisual.pPlot2D.sNumberTitle    = 'off';
            pVisual.pPlot2D.sMenuBar        = 'none';
            pVisual.pPlot2D.sName           = '2D';
            pVisual.pPlot2D.nOuterPosition  = [0 0 1 1];
            pVisual.pPlot2D.nPosition       = [0 0.04 1 0.935];
            pVisual.pPlot2D.sLegendLocation = 'NorthWest';
            pVisual.pPlot2D.sLegendTitle    = {'$\left|\Psi^m\right|$', 'Re $\Psi^m$'};
            pVisual.pPlot2D.sStyle          = sBaseStyle; % sStyle of plot: journal/print
            pVisual.pPlot2D.bSaveFigure     = bSaveFigure;
            pVisual.pPlot2D.sPath           = sSaveImagePath;
            pVisual.pPlot2D.sPrintType      = '-depsc2';
            pVisual.pPlot2D.sFilePrefix     = [ datestr(now, 'yyyymmdd_HHMMSS') '_' pTask.sDimension '_' pTask.sExampleName ];
            pVisual.pPlot2D.sFileName       = 'SOLUTION';
            switch pCalc.name
                case 'FEM'
                    pVisual.pPlot2D.sFileSuffix = [ pCalc.name '_N=' num2str(pCalc.N) ];
                case 'FFDS'
                    pVisual.pPlot2D.sFileSuffix = [ pCalc.name '_PARAM=' strrep(pVisual.thetaSymbolPrint,'/','-') ];
            end
            if strcmp(pCalc.type,'SSP')
                pVisual.pPlot2D.sFileSuffix = [ pCalc.type '_' pVisual.pPlot2D.sFileSuffix ];
            end
            pVisual.pPlot2D.sFileSuffix = [ pVisual.pPlot2D.sFileSuffix ',J=' num2str(pCalc.n-1) ',M=' num2str(pCalc.m-1) ];
            pVisual.pPlot2D.sFileExtension  = 'eps';
            pVisual.pPlot2D.bPlotPotential        = bPlotPotential;
            pVisual.pPlot2D.nDenominator          = max(abs(V(pTask, pVisual.pPlot2D.nPlotX)));
            pVisual.pPlot2D.sLegendTitlePotential = {'$\chi_I$'}; % {'$\frac{V}{\max\left|V\right|}$'}; % 
            pVisual.pPlot2D.sLineTypePotential    = '--c';
            pVisual.pPlot2D.nLineWidthPotential   = 2;
            if isfield(pTask, 'V_Q')
                pVisual.pPlot2D.sFileSuffix = [ 'Q=' num2str(int32(max(pTask.V_Q))) '_' pVisual.pPlot2D.sFileSuffix ];
            end
            switch pCalc.name
                case 'FEM'
                    pVisual.pPlot2D.sLabelY = [ pCalc.name ', n=' num2str(pCalc.N)];
                case 'FFDS'
                    pVisual.pPlot2D.sLabelY = [ pCalc.name ', '  '\theta=' thetaSymbolPrint];
            end
            if pVisual.bPlayMovie || pVisual.bSaveMovie
                pVisual.pPlot2D.pMovie.bResave      = false;
                pVisual.pPlot2D.pMovie.sName        = 'Solution';
                pVisual.pPlot2D.pMovie.sExtension   = 'avi';
                pVisual.pPlot2D.pMovie.sCompression = 'none';
                pVisual.pPlot2D.pMovie.nQuality     = 100;
                pVisual.pPlot2D.pMovie.nFPS         = 15;
                pVisual.pPlot2D.pMovie.sFullName    = [ pCalc.sWorkDirectory...
                    '\' pVisual.pPlot2D.pMovie.sName...
                    '.' pVisual.pPlot2D.pMovie.sExtension ];
                pVisual.pPlot2D.pMovie.bExist = ~~exist( pVisual.pPlot2D.pMovie.sFullName, 'file' );
                if pVisual.pPlot2D.pMovie.bExist
                    if pVisual.pPlot2D.pMovie.bResave
                        delete(pVisual.pPlot2D.pMovie.sFullName);
                        pVisual.pPlot2D.pMovie.bSave = true ;
                    else
                        pVisual.pPlot2D.pMovie.bSave = false;
                    end
                else
                    pVisual.pPlot2D.pMovie.bSave = pVisual.bSaveMovie;
                end
            end
            if pVisual.bPlotBehaviour3D
                pVisual.pPlot3D.sType = 'modulus'; % type of plot: modulus/real/imaginary
                switch pVisual.pPlot3D.sType
                    case 'modulus'
                        pVisual.pPlot3D.nLimitZ = [0 1];
                        pVisual.pPlot3D.sTitle  = '{\bf Modulus of ';
                        pVisual.pPlot3D.sLabelZ = '|\Psi(x,t)|';
                    case 'real'
                        pVisual.pPlot3D.sTitle  = '{\bf Real part of';
                        pVisual.pPlot3D.nLimitZ = [-1 1];
                        pVisual.pPlot3D.sLabelZ = 'Re(x,t)';
                    case 'imaginary'
                        pVisual.pPlot3D.sTitle  = '{\bf Imaginary part of';
                        pVisual.pPlot3D.nLimitZ = [-1 1];
                        pVisual.pPlot3D.sLabelZ = 'Im(x,t)';
                end
                pVisual.pPlot3D.sTitle = [pVisual.pPlot3D.sTitle ' solution of the 1D Schr\"{o}dinger equation.}' ];
                pVisual.pPlot3D.bPlayPause      = true;
                pVisual.pPlot3D.sLabelX         = 'x';
                pVisual.pPlot3D.sLabelY         = 't';
                pVisual.pPlot3D.nFontSize       = 14;
                pVisual.pPlot3D.sFontWeight     = 'bold';
                pVisual.pPlot3D.nLimitX         = [0 pCalc.X_0 ];
                pVisual.pPlot3D.nLimitY         = [0 pCalc.Tmax];
                pVisual.pPlot3D.nMeshGridX      = pCalc.x_mod;
                pVisual.pPlot3D.nMeshGridXI     = 0:(2*pCalc.h_mod):pCalc.X_0;
                pVisual.pPlot3D.nMeshGridY      = 0:pCalc.tau:pCalc.Tmax;
                pVisual.pPlot3D.nMeshGridYI     = 0:(2*pCalc.tau):pCalc.Tmax;
                pVisual.pPlot3D.nViewAngle      = [1 -4 1];
                pVisual.pPlot3D.nCameraPosition = [0 0 1];
                pVisual.pPlot3D.sUnits          = 'normalized';
                pVisual.pPlot3D.sNextPlot       = 'replacechildren';
                pVisual.pPlot3D.sNumberTitle    = 'off';
                pVisual.pPlot3D.sMenuBar        = 'none';
                pVisual.pPlot3D.sName           = '3D';
                pVisual.pPlot3D.nOuterPosition  = [0 0 1 1];
                pVisual.pPlot3D.nPosition       = [0 0.04 1 0.935];
                % to use 1/4 screen, issue:
                % pVisual.pPlot3D.nPosition = [ 1+nScreenSize(3)/2 nScreenSize(4)/2 nScreenSize(3)/2 nScreenSize(4)/2 ];
                % pVisual.pPlot3D.nCameraTarget    = [1 1 1];
                % pVisual.pPlot3D.nCameraUpVector  = [0 0 1];
                % pVisual.pPlot3D.nCameraViewAngle = [0 0 1];
                pVisual.pPlot3D.sStyle          = sBaseStyle;
                pVisual.pPlot3D.bSaveFigure     = bSaveFigure;
                pVisual.pPlot3D.sPath           = sSaveImagePath;
                pVisual.pPlot3D.sPrintType      = '-depsc2';
                pVisual.pPlot3D.sFilePrefix     = [ datestr(now, 'yyyymmdd_HHMMSS') '_' pTask.sDimension '_' pTask.sExampleName ];
                pVisual.pPlot3D.sFileName       = upper(pVisual.pPlot3D.sType);
                switch pCalc.name
                    case 'FEM'
                        pVisual.pPlot3D.sFileSuffix = [ 'N=' num2str(pCalc.N) ];
                    case 'FFDS'
                        pVisual.pPlot3D.sFileSuffix = ['PARAM=' strrep(pVisual.thetaSymbolPrint,'/','-')];
                end
                pVisual.pPlot3D.sFileExtension  = 'eps';
                pVisual.pPlot3D.bPlotPotential  = bPlotPotential;
                pVisual.pPlot3D.nDenominator    = max(abs(V(pTask, pVisual.pPlot3D.nMeshGridX)));
            end
        case '2D'
            pVisual.bPlotTimeBehaviour2D    = false;
            pVisual.pPlot3D.sType           = 'real'; % type of plot: modulus/real/imaginary
            switch pVisual.pPlot3D.sType
                case 'modulus'
                    pVisual.pPlot3D.nLimitZ = [0 1];
                    pVisual.pPlot3D.sTitle  = '{\bf Modulus of ';
                    pVisual.pPlot3D.sLabelZ = '|\Psi(x,y,t_0)|';
                case 'real'
                    pVisual.pPlot3D.sTitle  = '{\bf Real part of';
                    pVisual.pPlot3D.nLimitZ = [-1 1];
                    pVisual.pPlot3D.sLabelZ = 'Re(x,y,t_0)';
                case 'imaginary'
                    pVisual.pPlot3D.sTitle = '{\bf Imaginary part of';
                    pVisual.pPlot3D.nLimitZ = [-1 1];
                    pVisual.pPlot3D.sLabelZ = 'Im(x,y,t_0)';
            end
            pVisual.pPlot3D.sTitle = [pVisual.pPlot3D.sTitle ' solution of the 2D Schr\"{o}dinger equation.}' ];
            pVisual.pPlot3D.nPlotCount      = pVisual.nPlotCount;
            pVisual.pPlot3D.bPlayPause      = false;
            pVisual.pPlot3D.sLabelX         = 'x';
            pVisual.pPlot3D.sLabelY         = 'y';
            pVisual.pPlot3D.nFontSize       = 14;
            pVisual.pPlot3D.sFontWeight     = 'bold';
            pVisual.pPlot3D.nLimitX         = [0 pCalc.X_0];
            pVisual.pPlot3D.nLimitY         = [0 pTask.Y  ];
            pVisual.pPlot3D.nMeshGridX      = pCalc.x_mod;
            pVisual.pPlot3D.nMeshGridXI     = 0:(4*pCalc.h_mod):pCalc.X_0;
            pVisual.pPlot3D.nMeshGridY      = pCalc.delta:pCalc.delta:(pTask.Y-pCalc.delta);
            pVisual.pPlot3D.nMeshGridYI     = pCalc.delta:(4*pCalc.delta):(pTask.Y-pCalc.delta);
            pVisual.pPlot3D.nViewAngle      = [1 -4 1];
            pVisual.pPlot3D.nCameraPosition = [0 0 1];
            pVisual.pPlot3D.sUnits          = 'normalized';
            pVisual.pPlot3D.sNextPlot       = 'replacechildren';
            pVisual.pPlot3D.sNumberTitle    = 'off';
            pVisual.pPlot3D.sMenuBar        = 'none';
            pVisual.pPlot3D.sName           = '3D';
            pVisual.pPlot3D.nPosition       = [0 0.04 1 0.935];
            if pVisual.bPlayMovie || pVisual.bSaveMovie
                pVisual.pPlot3D.pMovie.bResave      = false;
                pVisual.pPlot3D.pMovie.sName        = 'Solution';
                pVisual.pPlot3D.pMovie.sExtension   = 'avi';
                pVisual.pPlot3D.pMovie.sCompression = 'none';
                pVisual.pPlot3D.pMovie.nQuality     = 100;
                pVisual.pPlot3D.pMovie.nFPS         = 15;
                pVisual.pPlot3D.pMovie.sFullName    = [ pCalc.sWorkDirectory...
                    '\' pVisual.pPlot3D.pMovie.sName...
                    '.' pVisual.pPlot3D.pMovie.sExtension ];
                pVisual.pPlot3D.pMovie.bExist = ~~exist( pVisual.pPlot3D.pMovie.sFullName, 'file' );
                if pVisual.pPlot3D.pMovie.bExist
                    if pVisual.pPlot3D.pMovie.bResave
                        delete(pVisual.pPlot3D.pMovie.sFullName);
                        pVisual.pPlot3D.pMovie.bSave = true ;
                    else
                        pVisual.pPlot3D.pMovie.bSave = false;
                    end
                else
                    pVisual.pPlot3D.pMovie.bSave = pVisual.bSaveMovie;
                end
            end
            pVisual.pPlot3D.sStyle          = sBaseStyle;
            pVisual.pPlot3D.bSaveFigure     = bSaveFigure;
            pVisual.pPlot3D.sPath           = sSaveImagePath;
            pVisual.pPlot3D.sPrintType      = '-depsc2';
            pVisual.pPlot3D.sFilePrefix     = [ datestr(now, 'yyyymmdd_HHMMSS') '_' pTask.sDimension '_' pTask.sExampleName ];
            pVisual.pPlot3D.sFileName       = upper(pVisual.pPlot3D.sType);
            switch pCalc.name
                case 'FEM'
                    pVisual.pPlot3D.sFileSuffix = [ 'nx=' num2str(pCalc.N) ',ny=' num2str(pCalc.eta) ];
                case 'FFDS'
                    pVisual.pPlot3D.sFileSuffix = [ 'theta=' strrep(pVisual.thetaSymbolPrint,'/','-') '_'...
                        'eta='   strrep(pVisual.etaSymbolPrint,  '/','-') ];
            end
            if strcmp(pCalc.type,'SSP')
                pVisual.pPlot3D.sFileSuffix = [ pCalc.type '_' pVisual.pPlot3D.sFileSuffix ];
            end
            pVisual.pPlot3D.sFileSuffix = [ pVisual.pPlot3D.sFileSuffix ',J=' num2str(pCalc.n-1) ',K=' num2str(pCalc.K-1) ',M=' num2str(pCalc.m-1) ];
            pVisual.pPlot3D.sFileExtension  = 'eps';
            pVisual.pPlot3D.bPlotPotential  = bPlotPotential;
            pVisual.pPlot3D.nDenominator    = max(abs(V(pTask, pVisual.pPlot3D.nMeshGridX)));
    end
else
    pVisual.bPlotTimeBehaviour2D = false;
    pVisual.nPlotCount = 0;
end
if pVisual.bPlotErrorTimeBehaviour || pVisual.bPlotErrorModulusTimeBehaviour
    pVisual.pPlotError.bResave    = false;
    pVisual.pPlotError.sName      = 'Error';
    pVisual.pPlotError.sExtension = 'eps';
    pVisual.pPlotError.sFullName  = [ pCalc.sWorkDirectory...
        '\' pVisual.pPlotError.sName...
        '.' pVisual.pPlotError.sExtension ];
    pVisual.pPlotError.bExist = ~~exist( pVisual.pPlotError.sFullName, 'file' );
    if pVisual.pPlotError.bExist
        if pVisual.pPlotError.bResave
            delete(pVisual.pPlotError.sFullName);
            pVisual.pPlotError.bSave = true ;
        else
            pVisual.pPlotError.bSave = false;
        end
    end
    pVisual.pPlotError.nFontSize       = 24;
    pVisual.pPlotError.sFontWeight     = 'bold';
    pVisual.pPlotError.sLineType       = {'--r', '-b'};
    pVisual.pPlotError.nLineWidth      = 3;
    pVisual.pPlotError.sTitle          = {'Absolute error', 'Relative error'};
    pVisual.pPlotError.sUnits          = 'normalized';
    pVisual.pPlotError.sNextPlot       = 'replacechildren';
    pVisual.pPlotError.sNumberTitle    = 'off';
    pVisual.pPlotError.sMenuBar        = 'none';
    pVisual.pPlotError.sName           = 'Error';
    pVisual.pPlotError.nPosition       = [0 0.04 1 0.935];
    pVisual.pPlotError.nLimitX         = [1 length(pCalc.nFreq)]; % [0 pCalc.m-1];
    pVisual.pPlotError.nFreq           = (pCalc.nFreq-1)*pCalc.tau;
    pVisual.pPlotError.nMaxValue       = pCalc.Tmax;
    pVisual.pPlotError.bAlignYLim      = false;
    pVisual.pPlotError.sLegendLocation = 'NorthWest';
    pVisual.pPlotError.sLegendTitle    = {'$C$', '$L_2$'};
    pVisual.pPlotError.sStyle          = sBaseStyle;
    pVisual.pPlotError.bSaveFigure     = bSaveFigure;
    pVisual.pPlotError.sPath           = sSaveImagePath;
    pVisual.pPlotError.sPrintType      = '-depsc2';
    pVisual.pPlotError.sFilePrefix     = [ datestr(now, 'yyyymmdd_HHMMSS') '_' pTask.sDimension '_' pTask.sExampleName ];
    pVisual.pPlotError.sFileName       = 'ERROR';
    pVisual.pPlotError.sFileExtension  = 'eps';
end
switch pTask.sDimension
    case '1D'
        switch pCalc.name
            case 'FEM'
                pVisual.pPlotError.sLabelY = [ pCalc.name ', n=' num2str(pCalc.N)];
                pVisual.pPlotError.sFileSuffix = [ 'N=' num2str(pCalc.N) ];
            case 'FFDS'
                pVisual.pPlotError.sLabelY = [ pCalc.name ', \theta=' thetaSymbolPrint];
                pVisual.pPlotError.sFileSuffix = [ 'PARAM=' strrep(pVisual.thetaSymbolPrint,'/','-')];
        end
        pVisual.pPlotError.sFileSuffix = [ pVisual.pPlotError.sFileSuffix ',J=' num2str(pCalc.n-1) ',M=' num2str(pCalc.m-1) ];
    case '2D'
        switch pCalc.name
            case 'FEM'
                pVisual.pPlotError.sLabelY = [ pCalc.name ', n=' num2str(pCalc.N) ', \eta=' etaSymbolPrint];
                pVisual.pPlotError.sFileSuffix = [ 'nx=' num2str(pCalc.N) ', ny=' num2str(pCalc.eta) ];
            case 'FFDS'
                pVisual.pPlotError.sLabelY = [ pCalc.name ', ' ...
                    '(\theta,\eta)=(' thetaSymbolPrint ', ' etaSymbolPrint ')' ];
                pVisual.pPlotError.sFileSuffix = [ 'PARAM1=' strrep(pVisual.thetaSymbolPrint,'/','-') ',PARAM2=' strrep(pVisual.etaSymbolPrint,'/','-') ];
        end
        pVisual.pPlotError.sFileSuffix = [ pVisual.pPlotError.sFileSuffix ',J=' num2str(pCalc.n-1) ',K=' num2str(pCalc.K-1) ',M=' num2str(pCalc.m-1) ];
end
if strcmp(pCalc.type,'SSP')
    pVisual.pPlotError.sFileSuffix = [ pCalc.type '_' pVisual.pPlotError.sFileSuffix ];
end
pVisual.pPlotError.sFileSuffix = [ pVisual.pPlotError.sFileSuffix ',FREQ=' num2str(nFrequency)];