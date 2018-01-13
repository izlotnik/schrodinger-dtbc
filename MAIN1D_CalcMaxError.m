clearContext; format short; echo off; mnct = maxNumCompThreads(1);

[pTask, pInitWave] = setTaskParams('Schrodinger', '1D', 'EX01r', getDataPath(), '');
% [pTask, pInitWave] = setTaskParams('Schrodinger', '1D', 'EX20r', getDataPath(), 'LAST');
pCalc = setCalcParams(pTask, {'FEM'}, [ 1 10 10 NaN 1 ], {'QR'}); T = pCalc.Tmax; clear pCalc;
% error = containers.Map({'calc', 'etalon', 'frequency', 'interval'}, {true, 'global', 1, [0 T; T T; 0 T/8; T/8 T; 0 T/4; T/4 T; 0 T/2; T/2 T]});
% error = containers.Map({'calc', 'etalon', 'frequency'}, {true, 'exact', 1}); % V=0
error = containers.Map({'calc'}, {false});
plotting = containers.Map({'solution', 'potential', 'count', 'solution_norm', 'norm', 'details'}, {false, false, 0, false, false, false});
solution = containers.Map('KeyType', 'char', 'ValueType', 'any');
% solution('type') = 'FFDS';
solution('type') = 'FEM';
switch pTask.sExampleName
    case {'EX01', 'EX01r'}
        nM = {300*([1 2 5 10])}; % {300*([5 10])}; % {300*[1 2]}; %  {300*([1 2 5 10])}; % {1500}; % {3000}; % {300*(1:6)}; %  {300*(1:5)}; %
        nR = {1:4};              % {[1 4]};        % {3};         %  {[1 4]};     %  
        switch solution('type')
            case 'FEM'
                nN = {5};             % {5};             % {5};    % {9};       % {9};         % {1:6};       % {3:9};       % {1:6};       % {9};         % {9};         % {3};   %  {5};        % {1:5};       % {1:10}
                nJ = {60*[1 2 5 10]}; % {60*[5 10]}; %  {5*60}; % {10*(2:9)};% {10*(2:15)}; % {50*(3:12)}; % {10*(2:15)}; % {90:30:780}; % {10*(1:6)};  % {30*(1:10)}; % {420}; % {60:30:780}; % {10:5:60};
                % ylim([4*10^(-5) 1.5*10^(0)]); % for absolute L2-norm
                % {'BestOutside', 'NorthEast'}
            case 'FFDS'
                nN = {[1/4 1/6 1/12 0]}; % {[1/4 -1/12 -1/6 -1/4 -1/2]}; % {[1/4 1/2 1 2]}; %
                nJ = {5000:1000:15000};  % {100:100:1000}; % {1000:500:5000}; % {1000:1000:5000}; %
        end
        % etalon = containers.Map({'type', 'N', 'J', 'M'}, {'FEM', 9, 90, 4*3000});
        etalon = containers.Map({'type', 'N', 'J', 'M', 'r'}, {'FEM', 9, 90, 4*3000, 4});
        % etalon = containers.Map({'type', 'N', 'J', 'M'}, {'FEM', 9, 60, 12*nM{1}(end)});
        % etalon = containers.Map({'type', 'N', 'J', 'M'}, {'FEM', 9, 30, 10*nM});
        etalon = prepareSolutionParams(etalon);
    case {'EX20', 'EX20r'}
        nM = {252*2.^(0:3)}; % {4*252}; % {16000}; % {16000*2.^(-6:0)};  % {500*2.^(0:5)};    %
        nR = {1:4};          % {3};     % {4};   % {4};     % {[1 2 4]};          % Richardson parameter
        switch solution('type')
            case 'FEM'
                nN = {9};       % {1:5};         % {6};       % {5};       % {1:2}
                nJ = {1*36};    % {(1:5)*36};    % {4*36};    % {4*36};    % {5*36};    %  (1:9:64)*5*36
            case 'FFDS'
                nN = {[1/4 1/6 1/12 0]};
                nJ = {(1:9:64)*5*36};
        end
        etalon = containers.Map({'type', 'N', 'J', 'M', 'r'}, {'FEM', 9, 4*36, 4*252*2^3, 4});
        % etalon = containers.Map({'type', 'N', 'J', 'M', 'r'}, {'FEM', 9, 4*nJ{1}(end), 4*nM{1}(end), 4});
        % etalon = containers.Map({'type', 'N', 'J', 'M', 'r'}, {'FEM', 9, 5*36, 4*32*252, 4});
        % etalon = containers.Map({'type', 'N', 'J', 'M', 'r'}, {'FEM', 9, 4*36, 4*16000, 4});
        % etalon = containers.Map({'type', 'N', 'J', 'M', 'r'}, {'FEM', 9, 2*36, 4*16000, 4});
        % etalon = containers.Map({'type', 'N', 'J', 'M'}, {'FEM', 5, 5*36, 12*16*1000});
        % etalon = containers.Map({'type', 'N', 'J', 'M'}, {'FEM', 5, 5*36, 12*nM});
        % etalon = containers.Map({'type', 'N', 'J', 'M'}, {'FEM', 5, 5*36, 10*nM});
        % etalon = containers.Map({'type', 'N', 'J', 'M'}, {'FEM', 9, 5*36, 50*nM});
        % etalon = containers.Map({'type', 'N', 'J', 'M'}, {'FEM', 9, 5*36, 1*nM});
        % etalon = containers.Map({'type', 'theta', 'J', 'M'}, {'FFDS', 1/12, 320*36, 10*nM});
        etalon = prepareSolutionParams(etalon);
    case {'EX22', 'EX22m'}
        nM = {4*4*32*36}; % {32*36*4*2.^(-5:0)}; %{32*36*4}; % {32*36}; % {1000}; % {4000}; %
        nR = {4};         % {1:4};        % {3}; % {3};       % {1};     % {1};    % {1};    %
        switch solution('type')
            case 'FEM'
                nN = {9};       % {9};       % {1:6};       % {5};        % {9};       % {1:5};
                nJ = {5*30};    % {2*30};    % {30*(1:10)}; % {30*(1:4)}; % {1*30};    % {30*(1:15)};
            case 'FFDS'
                nN = {[1/4 1/6 1/12 0]};    % {[1/4 1/2 1 2]};          % {[1/4 -1/12 -1/6 -1/4 -1/2]}; %
                nJ = {30*(1:15)};           % {linspace(450, 1350, 5)}; % {15*30*(1:9)}; %
        end
        % etalon = prepareSolutionParams(containers.Map({'type', 'N', 'J', 'M'}, {'FEM', 9, 450, 4*nM}));
        % etalon = prepareSolutionParams(containers.Map({'type', 'N', 'J', 'M', 'r'}, {'FEM', 9, 5*30, 4*4*32*36, 4}));
        etalon = prepareSolutionParams(containers.Map({'type', 'N', 'J', 'M', 'r'}, {'FEM', 9, 4*5*30, 4*4*32*36, 4}));
    case 'EX22r'
        nM = {4*72*2^5}; % {4*72*2.^(0:4)}; % {4*72*2.^(0:5)}; % 
        nR = {3};        % {1:4};           % {1:4};           % 
        switch solution('type')
            case 'FEM'
                nN = {1:6};       % {9};
                nJ = {30*(1:10)}; % {2*30};
        end
        etalon = prepareSolutionParams(containers.Map({'type', 'N', 'J', 'M', 'r'}, {'FEM', 9, 5*30, 4*4*32*72, 4}));
end
% error('solution') = etalon;

%% Calculate etalon solution
if error.isKey('etalon') && strcmp(error('etalon'), 'global')
    tic
    disp('Use global etalon solution...')
    error_e = containers.Map('KeyType', 'char', 'ValueType', 'any');
    error_e('calc') = false;
    error_e('frequency') =  GetValue(error, 'frequency', 1)*etalon('M')/max(nM{:});
    plotting_e = containers.Map({'solution', 'potential', 'count', 'solution_norm', 'norm', 'details'}, {false, false, 0, false, false, false});
    if etalon.isKey('r')
        fprintf('N=%d, J=%d, M=%d, r=%d\n', etalon('N'), etalon('J'), etalon('M'), etalon('r'));
        [ ~, ~, pVisual, pCalc, ~ ] = calcSolutionRichardson(pTask, pInitWave, etalon, error_e, plotting_e);
    else
        fprintf('N=%d, J=%d, M=%d\n', etalon('N'), etalon('J'), etalon('M'));
        [ ~, ~, pVisual, pCalc, ~ ] = calcSolutionCN(pTask, pInitWave, etalon, error_e, plotting_e);
    end
    if strcmp(etalon('type'), 'FFDS'), pCalc.N = 1; end
    error('U') = pCalc.U;
    error('params') = pCalc;
    % error('ylabel') = pVisual.pPlotError.sLabelY;
    error('suffix') = pVisual.pPlotError.sFileSuffix;
    toc
end

%%
for q=1:length(nJ)
    nMc = nM{q}; nRc = nR{q}; nJc = nJ{q}; nNc = nN{q};
    meca = NaN(length(nMc), length(nRc), length(nJc), length(nNc)); mela = meca; mecr = meca; melr = meca;
    eeca = meca; eela = meca; eecr = meca; eelr = meca;
    t = NaN(length(nMc), length(nRc), length(nJc), length(nNc), 11);
    iM = 1;
    s = 1; P = NaN(numel(meca), 4+11+8*size(GetValue(error, 'interval', [0 T; T T]), 1));
    tic
    for m=nMc
        iR = 1;
        for r=nRc
            iJ = 1;
            for n=nJc
                iN = 1;
                for N=nNc % Calc reference solution
                    solution('N') = N; solution('J') = n; solution('M') = m;
                    if ~ ( length(nRc) == 1 && r == 1 )
                        solution('r') = r;
                    end
                    solution = prepareSolutionParams(solution);
                    if solution.isKey('r')
                        [ rError, time, pVisual, pCalc, sFileName ] = calcSolutionRichardson(pTask, pInitWave, solution, error, plotting);
                    else
                        [ rError, time, pVisual, pCalc, sFileName ] = calcSolutionCN(pTask, pInitWave, solution, error, plotting);
                    end
                    % TIME
                    % t(iM, iR, iJ, iN) = time('total');
                    tt = [ time('total') time('solve') time('lhs') time('conv') time('rhs') time('others') ...
                        sum(time('kernel')) time('matrix') time('error') time('plot') time('file') ];
                    t(iM, iR, iJ, iN, :) = tt;
                    if ~isempty(rError)
                        xxx = parseErrors(rError, containers.Map({'T', 'M', 'P'}, {T, pCalc.m-1, GetValue(error, 'interval', [0 T; T T])}));
                        ela = subsref(xxx, struct('type', '()', 'subs', {'L2', {'abs'}, {'error'}}));
                        elr = subsref(xxx, struct('type', '()', 'subs', {'L2', {'rel'}, {'error'}}));
                        eca = subsref(xxx, struct('type', '()', 'subs', {'C',  {'abs'}, {'error'}}));
                        ecr = subsref(xxx, struct('type', '()', 'subs', {'C',  {'rel'}, {'error'}}));
                        tla = subsref(xxx, struct('type', '()', 'subs', {'L2', {'abs'}, {'time'}}));
                        tlr = subsref(xxx, struct('type', '()', 'subs', {'L2', {'rel'}, {'time'}}));
                        tca = subsref(xxx, struct('type', '()', 'subs', {'C',  {'abs'}, {'time'}}));
                        tcr = subsref(xxx, struct('type', '()', 'subs', {'C',  {'rel'}, {'time'}}));
                        % MAX
                        ind = 2; % 1;
                        mela(iM, iR, iJ, iN) = ela(ind); % mtla(iM, iR, iJ, iN) = tla(1);
                        meca(iM, iR, iJ, iN) = eca(ind); % mtca(iM, iR, iJ, iN) = tca(1);
                        melr(iM, iR, iJ, iN) = elr(ind); % mtlr(iM, iR, iJ, iN) = tlr(1);
                        mecr(iM, iR, iJ, iN) = ecr(ind); % mtcr(iM, iR, iJ, iN) = tcr(1);
                        % END
                        ind = 3; % 2; % 
                        eela(iM, iR, iJ, iN) = ela(ind); % etla(iM, iR, iJ, iN) = tla(2);
                        eeca(iM, iR, iJ, iN) = eca(ind); % etca(iM, iR, iJ, iN) = tca(2);
                        eelr(iM, iR, iJ, iN) = elr(ind); % etlr(iM, iR, iJ, iN) = tlr(2);
                        eecr(iM, iR, iJ, iN) = ecr(ind); % etcr(iM, iR, iJ, iN) = tcr(2);
                        %
                        fff = NaN(1, 8*size(GetValue(error, 'interval', [0 T; T T]), 1));
                        for u=1:size(GetValue(error, 'interval', [0 T; T T]), 1)
                            fff(((u-1)*8+1):(u*8)) = [ ela(u) eca(u) elr(u) ecr(u) tla(u) tca(u) tlr(u) tcr(u) ];
                        end
                        P(s, :) = [ N n m r tt fff ];
                    else
                        P(s, 1:15) = [ N n m r tt ];
                    end
                    s = s + 1;
                    iN = iN + 1;
                end
                iJ = iJ + 1;
            end
            iR = iR + 1;
        end
        iM = iM + 1;
    end
    toc
    % Plot errors
    if error.isKey('etalon')
        switch etalon('type')
            case 'FEM'
                sEtalonValue = [ 'N=' num2str(etalon('N')) ];
            case 'FFDS'
                sEtalonValue = [ 'PARAM=' strrep(sym2str(sym(etalon('theta'))),'/','-') ];
        end
        sFileSuffix = [ 'VS_' etalon('type') '_' sEtalonValue '_M=' num2str(etalon('M')) '_J=' num2str(etalon('J')) '_' ];
    else
        sFileSuffix = '';
    end
    pVisual.pPlotMaxError.sFileSuffix = [ sFileSuffix 'FREQ=' num2str(GetValue(error, 'frequency', 1)) ];
    for i={'abs', 'rel'} % {'abs'} %
        for j={'C', 'L_2'} % {'L_2'} %
            type = {i{1}, j{1}};
            if strcmp(type{1}, 'abs')
                if strcmp(type{2}, 'C'), e = meca; else e = mela; end
            else
                if strcmp(type{2}, 'C'), e = mecr; else e = melr; end
            end
            sFileName = plotMaxError(e, pVisual.pPlotMaxError, solution('type'), type, nMc, nRc, nJc, nNc, iM, iR, iJ, iN);
            csvwrite([ sFileName '.csv' ], e);
            sFileName = strrep(sFileName, [ '_' i{1} '_' j{1} '_' ], '_');
        end
    end
    %
    fprintf('Save: %s\n', [ sFileName '.csv' ]);
    csvwrite([ sFileName '.csv' ], P);
    save(sFileName, 'meca', 'mela', 'mecr', 'melr', 'eeca', 'eela', 'eecr', 'eelr', 't', ...
        'nM', 'nR', 'nJ', 'nN', 'nMc', 'nRc', 'nJc', 'nNc', 'iM', 'iR', 'iJ', 'iN', 'pVisual', 'solution');
end
[status, cmdout] = system('convertEPStoPDF.cmd'); maxNumCompThreads(mnct);
% P(:,[5 6 11 12])
% i=[1 3]; P(i+1,[5 6 11 12])./P(i,[5 6 11 12])-1