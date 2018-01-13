function [ pCalc, pFile ] = setCalcParams(pTask, sCalcName, Params, sCalcMethodName)
%% Solution parameters
bSaveSolution   = true ;
bLoadSolution   = true ;
bResaveSolution = false;
%% Error parameters
bCalcError      = true ; bCalcErrorModulus = false;
bSaveError      = true ;
bLoadError      = true ;
bResaveError    = false;
%%
switch pTask.name
    case 'Schrodinger'
        pCalc.name = sCalcName{1};
        switch pCalc.name
            case 'FFDS'
                pCalc.theta = Params(1);
            case 'FEM'
                pCalc.N = Params(1);
            otherwise
                logMessage('error', 'WrongCalcName', 'Unknown calculation name.\n Current value is %s', pCalc.name, mfilename);
        end
        pCalc.bStoreAllTimeLevels = Params(5);
        pCalc.method = sCalcMethodName{1};
        switch length(sCalcMethodName)
            case 1
                pCalc.type = 'CN'; % Crank-Nikolson (default value)
            case 2
                pCalc.type = sCalcMethodName{2};
            case 3
                pCalc.type = sCalcMethodName{2};
                pCalc.subtype = sCalcMethodName{3};
        end
        pCalc.isConstStep = true;               % TODO: add as input parameter
        pCalc.nBoundaryMultiplier = 1;          % TODO: add as input parameter
        pCalc.X_0 = pCalc.nBoundaryMultiplier * pTask.X;
        pCalc.pLBoundary.bUseDTBC = true;       % TODO: add as input parameter
        if pCalc.pLBoundary.bUseDTBC
            if length(sCalcName)==1
                pCalc.pLBoundary.sName = pCalc.name;
                switch pCalc.pLBoundary.sName
                    case 'FFDS'
                        pCalc.pLBoundary.theta = pCalc.theta;
                    case 'FEM'
                        pCalc.pLBoundary.N = pCalc.N;
                end
            else
                pCalc.pLBoundary.sName = sCalcName{2};
                switch pCalc.pLBoundary.sName
                    case 'FFDS'
                        pCalc.pLBoundary.theta = Params(4);
                    case 'FEM'
                        pCalc.pLBoundary.N = Params(4);
                end
            end
        end
        pCalc.pRBoundary.bUseDTBC = true; % TODO: add as input parameter
        if pCalc.pRBoundary.bUseDTBC
            if length(sCalcName)==1
                pCalc.pRBoundary.sName = pCalc.name;
                switch pCalc.pRBoundary.sName
                    case 'FFDS'
                        pCalc.pRBoundary.theta = pCalc.theta;
                    case 'FEM'
                        pCalc.pRBoundary.N = pCalc.N;
                end
            else
                pCalc.pRBoundary.sName = sCalcName{2};
                switch pCalc.pRBoundary.sName
                    case 'FFDS'
                        pCalc.pRBoundary.theta = Params(4);
                    case 'FEM'
                        pCalc.pRBoundary.N = Params(4);
                end
            end
        end
        pCalc.nInner = Params(2);
        pCalc.m      = Params(3)+1;
        switch pTask.sExampleName % TODO: add pCalc.Tmax as an input parameter
            case {'EX01','EX01a','EX01r'}
                pCalc.Tmax = 0.006;
            case {'EX02','EX02t'}
                pCalc.Tmax = 0.05;
            case 'EX02m'
                switch pTask.sDimension
                    case '1D'
                        pCalc.Tmax = 1*0.05;
                    case '2D'
                        pCalc.Tmax = 1*0.018; % 2011-CMMP: 0.02
                end
            case 'EX02n'
                pCalc.Tmax = 3*0.02;
            case 'EX02_VATRO'
                pCalc.Tmax = 0.25;
            case 'EX03'
                pCalc.Tmax = 0.03;    % 0.01;
            case 'EX04'
                pCalc.Tmax = 0.006;
            case 'EX05'
                pCalc.Tmax = 0.02;
            case 'EX10'
                pCalc.Tmax = 0.05;    % 0.02;
            case {'EX20', 'EX20r'}
                pCalc.Tmax = 16;
            case 'EX21'
                pCalc.Tmax = 0.05;
            case {'EX02k','EX22','EX22a','EX22k','EX22m','EX22n','EX23','EX31'}
                switch pTask.sDimension
                    case '1D'
                        pCalc.Tmax = 0.045;
                    case '2D'
                        pCalc.Tmax = 1*0.027; % EX22m: In order to plot mult by 5/3*
                end
            case {'EX22b', 'EX32'}
                pCalc.Tmax = 0.06;
            case {'EX22r', 'EX22s'}
                pCalc.Tmax = 2*0.045;
            case 'EX24'
                pCalc.Tmax = 2.5;
            case 'EX25'
                c_p = 1.054; c_m = 9.109;
                pCalc.Tmax = 4*90/(2*c_m/c_p);
            case 'EX30'
                pCalc.Tmax = 0.09;
            otherwise
                logMessage('error', 'WrongExampleName', 'Unknown example name.\n Current value is %s', pTask.sExampleName, mfilename);
        end
        pCalc.n = pCalc.nBoundaryMultiplier * pCalc.nInner + 1;
        switch pTask.sDimension % set values for dimentional parameters
            case '1D'
                if pCalc.isConstStep
                    pCalc.h   = pCalc.X_0 /(pCalc.n - 1);
                    pCalc.x   = 0:pCalc.h:pCalc.X_0;
                    pCalc.tau = pCalc.Tmax /(pCalc.m - 1);
                end
                if pCalc.X_0 < 0 || pCalc.Tmax < 0
                    logMessage('error', 'WrongDiscretizationValue', [ 'Values of X_0 and Tmax must be greather than zero.\n'...
                        ' Current values are %d and %d  respectively' ], [ pCalc.X_0 pCalc.Tmax], mfilename);
                end
                if pCalc.n < 2 || pCalc.m < 2
                    logMessage('error', 'WrongDiscretizationValue', [ 'Values of n and m must be greather than two.\n'...
                        ' Current values are %d and %d  respectively' ], [ pCalc.n pCalc.m], mfilename);
                end
            case '2D'
                pCalc.eta = Params(6);
                pCalc.K   = Params(7)+1;
                if pCalc.isConstStep
                    pCalc.h     = pCalc.X_0 /(pCalc.n - 1);
                    pCalc.delta = pTask.Y   /(pCalc.K - 1);
                    pCalc.x     = 0:pCalc.h:pCalc.X_0;
                    pCalc.tau   = pCalc.Tmax /(pCalc.m - 1);
                end
        end
        switch pCalc.name
            case 'FFDS'
                pCalc.n_mod         = pCalc.n;
                pCalc.h_mod         = pCalc.h;
                pCalc.n_per         = 1:pCalc.n;
                pCalc.x_mod         = pCalc.h*(pCalc.n_per-1);
                pCalc.n_perInner    = 1:(pCalc.nInner+1);
                pCalc.xInner        = pCalc.h_mod *(pCalc.n_perInner - 1);
            case 'FEM'
                pCalc.n_mod         =(pCalc.n-1)*pCalc.N+1;
                pCalc.h_mod         = pCalc.h / pCalc.N;
                pCalc.n_per         = 1:pCalc.n_mod;
                pCalc.x_mod         = 0:pCalc.h_mod:pCalc.X_0;
                pCalc.n_perInner    = 1:pCalc.N:(pCalc.nInner * pCalc.N + 1); % = 1:pCalc.N:pCalc.n_mod;
                pCalc.xInner        = pCalc.h_mod *(pCalc.n_perInner - 1);
            otherwise
                logMessage('error', 'WrongCalcName', 'Unknown calculation name.\n Current value is %s', pCalc.name, mfilename);
        end
    otherwise
        logMessage('error', 'UnknownTaskName', 'Unknown task name.\n Current value is %s', pTask.name, mfilename);
end
pCalc.nHashValue = calcHashValue(pCalc, pTask);
sBaseDirectory = [ pTask.sPath '\' pCalc.nHashValue ];
%% Create structure for info files with task and calc parameters
pFile.pInfo.sName      = 'Info';
pFile.pInfo.sExtension = 'txt';
pFile.pInfo.bSave      = true ;
pFile.pInfo.bResave    = false; % save if already exist
pFile.pInfo.bExist     = false;
pFile.pInfo.sFullName  = [ sBaseDirectory '\' pFile.pInfo.sName '.' pFile.pInfo.sExtension ];
sBaseRunID = pTask.sBaseRunID;
[ ~, ~, sMessageid ] = mkdir(sBaseDirectory);
if strcmp(sMessageid, 'MATLAB:MKDIR:DirectoryExists') % base directory already exist
    if strcmp(sBaseRunID, 'LAST') % get max run identifier
        y0 = dir(sBaseDirectory);
        y1 = y0(cellfun(@(x) x>0, {y0(:).isdir}) );
        y2 = y1(length({y1(:).name}));
        sBaseRunID = y2.name;
        if strcmp(sBaseRunID, '.') || strcmp(sBaseRunID, '..')
            sBaseRunID = '';
        end
    end
else
    sBaseRunID = '';
end
pFile.pInfo.bExist = ~~exist(pFile.pInfo.sFullName, 'file');
if pFile.pInfo.bExist
    if pFile.pInfo.bResave
        delete(pFile.pInfo.sFullName);
        pFile.pInfo.bSave = true;
    else
        pFile.pInfo.bSave = false;
    end
else
    pFile.pInfo.bSave = true;
end
if pFile.pInfo.bSave
    saveParameters(pFile.pInfo.sFullName, pTask);
    saveParameters(pFile.pInfo.sFullName, pCalc);
end
if isempty(sBaseRunID) % calc run identifier
    sRunID = datestr(now, 'yyyymmdd_HHMMSS_FFF');
else
    sRunID = sBaseRunID;
end
sWorkDirectory = [ sBaseDirectory '\' sRunID ];
pFile.pSolution.bExist    = false;
pFile.pSolution.bSave     = true ;
pFile.pSolution.sFullName = '';
if bSaveSolution || bLoadSolution % Create structure for file of solution
    pFile.pSolution.sName          = 'Solution';
    pFile.pSolution.sExtension     = 'mat';
    pFile.pSolution.bSave          = bSaveSolution;
    pFile.pSolution.bLoad          = bLoadSolution;
    pFile.pSolution.bResave        = bResaveSolution;
    pFile.pSolution.bExist         = false;
    pFile.pSolution.sFullName      = [ sWorkDirectory '\' pFile.pSolution.sName '.' pFile.pSolution.sExtension ];
end
pFile.pError.bSave            = bSaveError;
pFile.pError.sFullName        = '';
pFile.pErrorModulus.bSave     = bSaveError;
pFile.pErrorModulus.sFullName = '';
if bSaveError || bLoadError % Create structure for error's file
    pFile.pError.sName             = 'Error';
    pFile.pError.sExtension        = 'mat';
    pFile.pError.bLoad             = bLoadError;
    pFile.pError.bResave           = bResaveError;
    pFile.pError.bExist            = false;
    pFile.pError.sFullName         = [ sWorkDirectory '\' pFile.pError.sName  '.' pFile.pError.sExtension ];
    %
    pFile.pErrorModulus.sName      = 'ErrorModulus';
    pFile.pErrorModulus.sExtension = 'mat';
    pFile.pErrorModulus.bLoad      = bLoadError;
    pFile.pErrorModulus.bResave    = bResaveError;
    pFile.pErrorModulus.bExist     = false;
    pFile.pErrorModulus.sFullName  = [ sWorkDirectory '\' pFile.pErrorModulus.sName '.' pFile.pErrorModulus.sExtension ];
end
[ ~, ~, sMessageid ] = mkdir(sWorkDirectory);
if strcmp(sMessageid, 'MATLAB:MKDIR:DirectoryExists')
    if bSaveSolution || bLoadSolution
        pFile.pSolution.bExist = ~~exist(pFile.pSolution.sFullName, 'file');
        if pFile.pSolution.bExist
            if pFile.pSolution.bResave
                delete(pFile.pSolution.sFullName);
                pFile.pSolution.bSave = true;
                pFile.pSolution.bLoad = false;
            else
                pFile.pSolution.bSave = false;
            end
        else
            pFile.pSolution.bLoad = false;
        end
    end
    if bSaveError || bLoadError
        pFile.pError.bExist = ~~exist(pFile.pError.sFullName, 'file');
        if pFile.pError.bExist
            if pFile.pError.bResave
                delete(pFile.pError.sFullName);
                pFile.pError.bSave = true;
                pFile.pError.bLoad = false;
            else
                pFile.pError.bSave = false;
            end
        else
            pFile.pError.bLoad = false;
        end
        pFile.pErrorModulus.bExist = ~~exist(pFile.pErrorModulus.sFullName, 'file');
        if pFile.pErrorModulus.bExist
            if pFile.pErrorModulus.bResave
                delete(pFile.pErrorModulus.sFullName);
                pFile.pErrorModulus.bSave = true;
                pFile.pErrorModulus.bLoad = false;
            else
                pFile.pErrorModulus.bSave = false;
            end
        else
            pFile.pErrorModulus.bLoad = false;
        end
    end
end
pCalc.bCalcError = bCalcError;
pCalc.bCalcErrorModulus = bCalcErrorModulus;
%
pCalc.sBaseDirectory = sBaseDirectory;
pCalc.sWorkDirectory = sWorkDirectory;