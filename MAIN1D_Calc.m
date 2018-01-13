clearContext; format short; echo off;

% [pTask, pInitWave] = setTaskParams('Schrodinger', '1D', 'EX01r', getDataPath(), '');
[pTask, pInitWave] = setTaskParams('Schrodinger', '1D', 'EX22b', getDataPath(), 'LAST');
% [pTask, pInitWave] = setTaskParams('Schrodinger', '1D', 'EX20r', getDataPath(), 'LAST');
pCalc = setCalcParams(pTask, {'FEM'}, [ 1 10 10 NaN 1 ], {'QR'}); T = pCalc.Tmax; clear pCalc;
% error = containers.Map({'calc', 'etalon', 'frequency', 'interval'}, {true, 'global', 1, [0 T; T T; 0 T/8; T/8 T; 0 T/4; T/4 T; 0 T/2; T/2 T]});
% error = containers.Map({'calc', 'etalon', 'frequency'}, {true, 'exact', 1}); % V=0
error = containers.Map({'calc'}, {false}); % plot solution
plotting = containers.Map({'solution', 'potential', 'count', 'solution_norm', 'norm', 'details'}, {true, true, 0, false, false, false}); % plot solution
% plotting = containers.Map({'solution', 'potential', 'count', 'solution_norm', 'norm', 'details'}, {true, false, 0, false, false, false}); % V=0
% plotting = containers.Map({'solution', 'potential', 'count', 'solution_norm', 'norm', 'details'}, {false, false, 0, false, true, false});
% plotting = containers.Map({'solution', 'potential', 'count', 'solution_norm', 'norm', 'details'}, {false, false, 0, false, false, false});
switch pTask.sExampleName
    case {'EX01', 'EX01a', 'EX01r'}
        % solution = containers.Map({'type', 'theta', 'J', 'M, 'r'}, {'FFDS', 1/12, 100, 300, 4}); plotting('count') = 30;
        % solution = containers.Map({'type', 'theta', 'J', 'M'}, {'FFDS', 1/12, 1000, 3000});
        % solution = containers.Map({'type', 'N', 'J', 'M'}, {'FEM', 1, 5000, 3000});
        % solution = containers.Map({'type', 'N', 'J', 'M'}, {'FEM', 2, 300, 3000});
        % solution = containers.Map({'type', 'N', 'J', 'M', 'r'}, {'FEM', 5, 300, 300, 4});
        % solution = containers.Map({'type', 'N', 'J', 'M', 'r'}, {'FEM', 5, 1*30, 1*3000, 2});
        % solution = containers.Map({'type', 'N', 'J', 'M', 'r'}, {'FEM', 9, 30, 300, 4}); % plotting('count') = 30;
        % solution = containers.Map({'type', 'N', 'J', 'M'}, {'FEM', 9, 30, 3000});
        solution = containers.Map({'type', 'N', 'J', 'M'}, {'FEM', 5, 120, 3000});
        % solution = containers.Map({'type', 'N', 'J', 'M', 'r'}, {'FEM', 5, 1*30, 50*300, 4}); %% 2014/04/13
        % solution = containers.Map({'type', 'N', 'J', 'M'}, {'FEM', 9, 100, 3000});
        % solution = containers.Map({'type', 'N', 'J', 'M'}, {'FEM', 9, 10, 12*16*16});
        solution = prepareSolutionParams(solution);
        etalon = containers.Map({'type', 'N', 'J', 'M', 'r'}, {'FEM', 9, 30, 4*300, 4}); 
        % etalon = containers.Map({'type', 'N', 'J', 'M', 'r'}, {solution('type'), solution('N'), solution('J'), 4*solution('M'), 4}); 
        % etalon = containers.Map({'type', 'N', 'J', 'M'}, {'FEM', 9, 60, 12*solution('M')}); 
        etalon = prepareSolutionParams(etalon);
    case 'EX02m'
        solution = prepareSolutionParams(containers.Map({'type', 'N', 'J', 'M'}, {'FEM', 9, (7/12)*1200, 1000}));
        etalon   = prepareSolutionParams(containers.Map({'type', 'N', 'J', 'M'}, {'FEM', 9, 60, 10*solution('M')}));
    case {'EX20', 'EX20r'}
        % solution = containers.Map({'type', 'N', 'J', 'M'}, {'FEM', 9, 1*36, 8*252}); plotting('count') = 4*solution('M');
        solution = containers.Map({'type', 'N', 'J', 'M', 'r'}, {'FEM', 9, 1*36, 2*252, 3});
        plotting('count') = solution('r')*solution('M');
        % solution = containers.Map({'type', 'N', 'J', 'M', 'r'}, {'FEM', 9, 2*36, 4*252*4, 4});
        % solution = containers.Map({'type', 'N', 'J', 'M', 'r'}, {'FEM', 9, 2*36, 4*252*2, 3});
        % solution = containers.Map({'type', 'N', 'J', 'M', 'r'}, {'FEM', 9, 4*36, 4*252*2^3, 4});
        % solution = containers.Map({'type', 'N', 'J', 'M'}, {'FEM', 9, 1*36, 1*16000});
        % solution = containers.Map({'type', 'N', 'J', 'M'}, {'FEM', 9, 5*36, 16000});
        % solution = containers.Map({'type', 'N', 'J', 'M'}, {'FEM', 9, 160*36, 500});
        % solution = containers.Map({'type', 'N', 'J', 'M'}, {'FEM', 9, 160*36, 10*1600});
        solution = prepareSolutionParams(solution);
        etalon = containers.Map({'type', 'N', 'J', 'M', 'r'}, {'FEM', 9, 4*36, 4*252*2^3, 4});
        % etalon = containers.Map({'type', 'N', 'J', 'M', 'r'}, {'FEM', 9, 4*36, 4*4*252*2^3, 4});
        % etalon = containers.Map({'type', 'N', 'J', 'M', 'r'}, {'FEM', 9, 4*4*36, 4*252*2^3, 4});
        % etalon = containers.Map({'type', 'N', 'J', 'M', 'r'}, {'FEM', 9, 5*36, 4*32*252, 4});
        % etalon = containers.Map({'type', 'N', 'J', 'M', 'r'}, {'FEM', 9, 2*36, 4*16000, 4}); 
        % etalon = containers.Map({'type', 'N', 'J', 'M'}, {'FEM', 5, 320*36, 12*10*1600}); 
        % etalon = containers.Map({'type', 'theta', 'J', 'M'}, {'FFDS', 1/12, 320*36, 12*solution('M')});
        etalon = prepareSolutionParams(etalon);
    case {'EX22', 'EX22a', 'EX22m'}
        solution = prepareSolutionParams(containers.Map({'type', 'N', 'J', 'M', 'r'}, {'FEM', 9, 2*30, 32*36, 2})); 
        % solution = prepareSolutionParams(containers.Map({'type', 'N', 'J', 'M', 'r'}, {'FEM', 9, 5*30, 4*4*32*36, 4}));
        % solution = prepareSolutionParams(containers.Map({'type', 'N', 'J', 'M', 'r'}, {'FEM', 9, 2*30, 2*4*36, 4})); plotting('count') = 24;
        % solution = prepareSolutionParams(containers.Map({'type', 'N', 'J', 'M', 'r'}, {'FEM', 9, 2*30, 8*36, 3}));
        % solution = prepareSolutionParams(containers.Map({'type', 'theta', 'J', 'M'}, {'FFDS', 1/12, 1200, 1000}));
        % etalon   = prepareSolutionParams(containers.Map({'type', 'N', 'J', 'M'}, {'FEM', 9, 15*solution('J'), 4*solution('M')}));
        % etalon   = prepareSolutionParams(containers.Map({'type', 'N', 'J', 'M', 'r'}, {'FEM', 9, 5*solution('J'), 12*solution('M'), 4}));
        etalon   = prepareSolutionParams(containers.Map({'type', 'N', 'J', 'M', 'r'}, {'FEM', 9, 5*30, 4*4*32*36, 4}));
        % etalon = prepareSolutionParams(containers.Map({'type', 'N', 'J', 'M', 'r'}, {'FEM', 9, 4*5*30, 4*4*32*36, 4}));
        % etalon = prepareSolutionParams(containers.Map({'type', 'N', 'J', 'M', 'r'}, {'FEM', 9, 5*30, 4*4*4*32*36, 4}));
    case 'EX22b' % 2014/11/10
        solution = prepareSolutionParams(containers.Map({'type', 'N', 'J', 'M', 'r'}, {'FEM', 9, 300, 12000, 1}));
        plotting('count') = 600; % solution('M'); % 4*solution('M'); % 36; % 
        etalon   = prepareSolutionParams(containers.Map({'type', 'N', 'J', 'M', 'r'}, {'FEM', 9, 4*675, 4*600, 1}));
    case 'EX32' % 2014/11/10
        solution = prepareSolutionParams(containers.Map({'type', 'theta', 'J', 'M', 'r'}, {'FFDS', 1/12, 675, 600, 1}));
        plotting('count') = solution('M'); % 10; % 4*solution('M'); % 36; % 
        etalon   = prepareSolutionParams(containers.Map({'type', 'theta', 'J', 'M', 'r'}, {'FFDS', 1/12, 4*675, 4*600, 1}));    
    case 'EX22r'
        % solution = prepareSolutionParams(containers.Map({'type', 'N', 'J', 'M', 'r'}, {'FEM', 9, 2*30, 4*72*2^3, 2})); 
        % solution = prepareSolutionParams(containers.Map({'type', 'N', 'J', 'M', 'r'}, {'FEM', 9, 2*30, 4*72*2^5, 4})); 
        solution = prepareSolutionParams(containers.Map({'type', 'N', 'J', 'M', 'r'}, {'FEM', 9, 2*30, 8*72, 4})); 
        plotting('count') = 4*solution('M'); % 36; % 
        etalon   = prepareSolutionParams(containers.Map({'type', 'N', 'J', 'M', 'r'}, {'FEM', 9, 5*30, 4*4*32*72, 4}));
    case 'EX22s' % 2013/11/08
        % nM_F = 2.^(-3:1); % 1; %
        % nJ_F = 2.^( 1:6); % 1; %
        solution = prepareSolutionParams(containers.Map({'type', 'N', 'J', 'M'}, {'FEM', 9, 30, 4000}));
        etalon   = prepareSolutionParams(containers.Map({'type', 'N', 'J', 'M'}, {'FEM', 9, 4*solution('J'), 4*solution('M')}));
    case 'EX24' % Moxley
        % solution = prepareSolutionParams(containers.Map({'type', 'N', 'J', 'M'}, {'FEM', 9, 100*16, 50000})); % Moxley
        solution = prepareSolutionParams(containers.Map({'type', 'N', 'J', 'M'}, {'FEM', 9, 10*16, 2*50000}));
        etalon   = prepareSolutionParams(containers.Map({'type', 'N', 'J', 'M'}, {'FEM', 9, 10*16, 4*50000}));
    case 'EX25' % potential step by Sullivan - 2012 % 2013/11/08
        % nM_F = 2.^(-4:0); % 1; %
        % nJ_F = 2.^( 0:6); % 1; %
        solution = prepareSolutionParams(containers.Map({'type', 'N', 'J', 'M'}, {'FEM', 9, 20, 4*1000}));
        etalon   = prepareSolutionParams(containers.Map({'type', 'N', 'J', 'M'}, {'FEM', 9, 2*4*solution('J'), 2*4*solution('M')}));
end
% error('solution') = etalon;

%% Task properties
if isfield(pTask, 'V_x') % check values on the artificial boundaries
    fprintf('V(0)-V_L=%d, V(X)-V_R=%d\n', V(pTask, 0) - pTask.V_Linf, V(pTask, pTask.X) - pTask.V_Rinf);
    fprintf('max_V(abs(psi0))=%g\n', max(abs(psi0(pTask, pInitWave, pTask.V_x(1):pTask.V_x(end)))));
end

%% Define execution parameters
nM_F = 1; % 2.^(-3:1); % 2.^(-2:1); %  2.^(-4:0); % 2.^0; % 2.^(-2:1)
nJ_F = 1; % 2.^( 1:6); % 1; % 2.^( 0:6); % 2^4; % [1:0.5:6 6:3:60]; % [1:10 10:5:50 50:10:150]; % 2.^0; % [1:2 4:5]; % 2.^(-5:0) % 2.^(-4:0) % 1:10 %
lM = length(nM_F); lJ = length(nJ_F); P = NaN(lM*lJ, 4+11+8*size(GetValue(error, 'interval', [0 T; T T]), 1));
for m=1:lM, for j=1:lJ, P((m-1)*lJ+j, 1:2) = [ nJ_F(j)*solution('J') nM_F(m)*solution('M') ]; end, end

%% Calculate etalon solution
if error.isKey('etalon') && strcmp(error('etalon'), 'global')
    tic
    disp('Use global etalon solution...')
    error_e = containers.Map('KeyType', 'char', 'ValueType', 'any'); 
    error_e('calc') = false; 
    error_e('frequency') = error('frequency')*etalon('M')/max(P(:,2));
    plotting_e = containers.Map({'solution', 'potential', 'count', 'solution_norm', 'norm', 'details'}, {false, false, 0, false, false, false});
    if etalon.isKey('r')
        [ ~, ~, pVisual, pCalc, ~ ] = calcSolutionRichardson(pTask, pInitWave, etalon, error_e, plotting_e);
    else
        [ ~, ~, pVisual, pCalc, ~ ] = calcSolutionCN(pTask, pInitWave, etalon, error_e, plotting_e);
    end
    if strcmp(etalon('type'), 'FFDS'), pCalc.N = 1; end
    error('U') = pCalc.U;
    error('params') = pCalc;
    % error('ylabel') = pVisual.pPlotError.sLabelY;
    error('suffix') = pVisual.pPlotError.sFileSuffix;
    %
    fprintf('EPS0=max_{Omega/(0,X)} U(x,0)=%d\nEPS1=max_{[0,T]} U(0,t)=%d\nEPS2=max_{[0,T]} U(X,t)=%d\n', ...
        max(abs(pCalc.U([1 end], 1))), max(abs(pCalc.U(1, :))), max(abs(pCalc.U(end, :))));
    toc
else % display solution (initial wave) parameters
    fprintf('EPS0=max_{Omega/(0,X)} psi_G(x,0)=%d\n', max(abs(GaussianWave(pTask, pInitWave, 0, [0 pTask.X]))));
%     pp = 1; pCalc = setCalcParams(pTask, {solution('type')}, [solution('N') P(pp, 1) P(pp, 2) NaN 1], {solution('matrix')});
%     fprintf('EPS0=max_{Omega/(0,X)} psi_G(x,0)=%d\nEPS1=max_{[0,T]} psi_G(0,t)=%d\nEPS2=max_{[0,T]} psi_G(X,t)=%d\n', ...
%         max(abs(GaussianWave(pTask, pInitWave, 0, [0 pCalc.X_0]))), ...
%         max(abs(GaussianWave(pTask, pInitWave, pCalc.tau*((1:pCalc.m)-1), 0        ))), ...
%         max(abs(GaussianWave(pTask, pInitWave, pCalc.tau*((1:pCalc.m)-1), pCalc.X_0))));
end

%% Calculate solution
for pp=1:(lM*lJ)
    fprintf('pp=%d (%.2f)\n', pp, pp/(lM*lJ)*100);
    solution('J') = P(pp, 1);
    solution('M') = P(pp, 2);
    if error.isKey('etalon') && strcmp(error('etalon'), 'local')
        disp('Use local etalon solution...')
        % modify: for example change DTBC to SDTBC (use DTBC here and SDTBC below)
        solution_e = solution;
        remove(solution_e, 'r'); solution_e = prepareSolutionParams(solution_e);
        error_e    = containers.Map('KeyType', 'char', 'ValueType', 'any'); error_e('calc') = false; % error_e = error;
        plotting_e = containers.Map({'solution', 'potential', 'count', 'solution_norm', 'norm', 'details'}, {false, false, 0, false, false, false});
        if solution_e.isKey('r')
            [ ~, ~, pVisual, pCalc, ~ ] = calcSolutionRichardson(pTask, pInitWave, solution_e, error_e, plotting_e);
        else
            [ ~, ~, pVisual, pCalc, ~ ] = calcSolutionCN(pTask, pInitWave, solution_e, error_e, plotting_e);
        end
        if strcmp(solution_e('type'), 'FFDS'), pCalc.N = 1; end
        error('U') = pCalc.U;
        error('params') = pCalc;
        % error('ylabel') = pVisual.pPlotError.sLabelY;
        error('suffix') = 'DTBC_VS_SDTBC';
        %
        fprintf('EPS0=max_{Omega/(0,X)} U(x,0)=%d\nEPS1=max_{[0,T]} U(0,t)=%d\nEPS2=max_{[0,T]} U(X,t)=%d\n', ...
            max(abs(pCalc.U([1 end], 1))), max(abs(pCalc.U(1, :))), max(abs(pCalc.U(end, :))));
    end
    if solution.isKey('r')
        [ rError, time, pVisual, pCalc, sFileName ] = calcSolutionRichardson(pTask, pInitWave, solution, error, plotting);
    else
        [ rError, time, pVisual, pCalc, sFileName ] = calcSolutionCN(pTask, pInitWave, solution, error, plotting);
    end
    % TIME
    tt = [ time('total') time('solve') time('lhs') time('conv') time('rhs') time('others') ...
        sum(time('kernel')) time('matrix') time('error') time('plot') time('file') ];
    if solution.isKey('r')
        r = solution('r');
    else
        r = 1;
    end
    switch solution('type')
        case 'FEM'
            param = solution('N');
        case 'FFDS'
            param = solution('theta');
    end
    if ~isempty(rError)
%         % MAX
%         [mla, ila] = max(rError.L2Abs); [mca, ica] = max(rError.CAbs);
%         [mlr, ilr] = max(rError.L2Rel); [mcr, icr] = max(rError.CRel);
%         % END
%         ela = rError.L2Abs(end); eca = rError.CAbs(end);
%         elr = rError.L2Rel(end); ecr = rError.CRel(end);
%         %
%         % subsref(x, struct('type', '()', 'subs', {'L2', {'abs'}, {'time'}}))
%         subsref(xxx, struct('type', '()', 'subs', {'L2', {'abs'}, {'error'}}))-[mla; ela]
%         subsref(xxx, struct('type', '()', 'subs', {'L2', {'rel'}, {'error'}}))-[mlr; elr]
%         subsref(xxx, struct('type', '()', 'subs', {'C', {'abs'}, {'error'}}))-[mca; eca]
%         subsref(xxx, struct('type', '()', 'subs', {'C', {'rel'}, {'error'}}))-[mcr; ecr]
%         subsref(xxx, struct('type', '()', 'subs', {'L2', {'abs'}, {'time'}}))-[(ila-1)/(pCalc.m-1)*pCalc.Tmax; pCalc.Tmax]
%         subsref(xxx, struct('type', '()', 'subs', {'L2', {'rel'}, {'time'}}))-[(ilr-1)/(pCalc.m-1)*pCalc.Tmax; pCalc.Tmax]
%         subsref(xxx, struct('type', '()', 'subs', {'C', {'abs'}, {'time'}}))-[(ica-1)/(pCalc.m-1)*pCalc.Tmax; pCalc.Tmax]
%         subsref(xxx, struct('type', '()', 'subs', {'C', {'rel'}, {'time'}}))-[(icr-1)/(pCalc.m-1)*pCalc.Tmax; pCalc.Tmax] 
        xxx = parseErrors(rError, containers.Map({'T', 'M', 'P'}, {T, pCalc.m-1, GetValue(error, 'interval', [0 T; T T])}));
        ela = subsref(xxx, struct('type', '()', 'subs', {'L2', {'abs'}, {'error'}}));
        elr = subsref(xxx, struct('type', '()', 'subs', {'L2', {'rel'}, {'error'}}));
        eca = subsref(xxx, struct('type', '()', 'subs', {'C',  {'abs'}, {'error'}}));
        ecr = subsref(xxx, struct('type', '()', 'subs', {'C',  {'rel'}, {'error'}}));
        tla = subsref(xxx, struct('type', '()', 'subs', {'L2', {'abs'}, {'time'}}));
        tlr = subsref(xxx, struct('type', '()', 'subs', {'L2', {'rel'}, {'time'}}));
        tca = subsref(xxx, struct('type', '()', 'subs', {'C',  {'abs'}, {'time'}}));
        tcr = subsref(xxx, struct('type', '()', 'subs', {'C',  {'rel'}, {'time'}}));
        %
        fff = NaN(1, 8*size(GetValue(error, 'interval', [0 T; T T]), 1));
        for u=1:size(GetValue(error, 'interval', [0 T; T T]), 1)
            fff(((u-1)*8+1):(u*8)) = [ ela(u) eca(u) elr(u) ecr(u) tla(u) tca(u) tlr(u) tcr(u) ];
        end
        P(pp, :) = [ param P(pp, 1) P(pp, 2) r tt fff ];
    else
       P(pp, 1:15) = [ param P(pp, 1) P(pp, 2) r tt ];
    end
end
fprintf('Save: %s\n', [ sFileName '.csv' ]);
csvwrite([ sFileName '.csv' ], P);

%%
if pVisual.bPlotBehaviour3D % Plot and save plot of all time-levels solution
    pVisual.pPlot3D = initGraphics('3D', pVisual.pPlot3D);
    plot3DSolutionBehaviour(transpose(pCalc.U), pVisual.pPlot3D, pTask); set(gcf, 'Renderer', 'zbuffer');
    saveImage(gcf, pVisual.pPlot3D, '', mfilename); close(gcf);
end
if pVisual.bPlayMovie % Play movie
    initGraphics('2D', pVisual.pPlot2D);
    movie(pVisual.pPlot2D.pMovie.pFrame);
    % for mm=1:pCalc.m, image(pVisual.pPlot2D.pMovie.pFrame(mm).cdata), colormap(pVisual.pPlot2D.pMovie.pFrame(mm).colormap), pause, end
end
if pVisual.pPlot2D.pMovie.bSave % Save movie
    saveMovie(pVisual.pPlot2D.pMovie.pFrame, pVisual.pPlot2D.pMovie);
    % nQuality = 100; sCodec =  'yv12'; %'divx'; % 'i420'; % 'Indeo5'; % 'xvid'; % 'iyuv'; % 'none'; %
    % movie2avi(pVisual.pPlot2D.pMovie.pFrame, ...
    %    [ pVisual.pPlot2D.pMovie.sFullName '_' sCodec '_' num2str(nQuality) '.avi' ], 'compression', sCodec, 'quality', nQuality);
end

%% Plot solution behaviour
if false % true % 
    [pVisual, pCalc] = setVisualParams(pTask, pCalc, ...
        [1 plotting('potential') plotting('count') 0 0], ...
        [1 GetValue(error, 'frequency', 1) 1 0], [0 1]);
    % pVisual.pPlot2D = rmfield(pVisual.pPlot2D, 'nLimitY');
    % EX22r:
    pVisual.pPlot2D.nLimitY = [-2.15 2.15];
    % EX24:
    % pVisual.pPlot2D.nLineWidth = 1;
    % pVisual.pPlot2D.nLimitY = [-0.25 0.25];
    if pVisual.bPlotTimeBehaviour % Plot solution behaviour
        pVisual.pPlot2D = initGraphics('2D', pVisual.pPlot2D);
        %
        params = error;
        if params.isKey('etalon') && ( strcmp(params('etalon'), 'global') || strcmp(params('etalon'), 'local') )
            useEtalon = true;
            U_E = params('U');
            pCalc_E = params('params');
        else
            useEtalon = false;
        end
        switch pCalc.name
            case 'FFDS'
                nn = pCalc.n_perInner; xx = pCalc.xInner;
            case 'FEM'
                nn = pCalc.n_per;      xx = pCalc.x_mod;
        end
        if useEtalon
            s = (size(U_E, 2)-1)/(pCalc.m-1); if s~=int32(s), keyboard, end
            x = (pCalc.n_per-1)'*pCalc.h_mod;
        end
        for k=1:pCalc.m % 5 % 
            if mod(k-1, (pCalc.m - 1) / pVisual.pPlot2D.nPlotCount) == 0 && pVisual.bPlotTimeBehaviour
                U = pCalc.U(nn, k);
                if useEtalon
                    if pCalc_E.n == pCalc.n && pCalc_E.N == 1 && ~isfield(pCalc, 'N')
                        y = U_E(:, (k-1)*s+1);
                    else
                        D = PFD(U_E(:, (k-1)*s+1), pCalc_E.n-1, pCalc_E.N);
                        y = PN(D, pCalc_E.N, pCalc_E.h, x, pCalc_E.n-1);
                    end
                else
                    y = transpose(GaussianWave(pTask, pInitWave, pCalc.tau*(k-1), xx));
                end
                if isfield(pTask, 'V_x')
                    nV = nn( xx >= pTask.V_x(1) & xx <= pTask.V_x(end) );
                    fprintf('m=%d: max_V(abs(U))=%g, max_V(abs(E))=%g, max_V(abs(U-E))=%g\n', k-1, ...
                        max(abs(U(nV, 1))), max(abs(y(nV, 1))), max(abs(U(nV, 1)-y(nV, 1))));
                end
                % U = U - y; % plot error
                pVisual.pPlot2D = plot2DSolutionBehaviour(U, pVisual.pPlot2D, pTask, k);
                saveImage(gcf, pVisual.pPlot2D, ['M=' num2str(k-1)], mfilename);
            end
        end
    end
end
[status, cmdout] = system('convertEPStoPDF.cmd');