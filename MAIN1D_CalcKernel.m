clearContext; format short; echo off;
%
sExampleName   = 'EX01';    %'EX25';    % 'EX23'; %
sDataDirectory = 'D:\Matlab\R2011b';   % 'E:';      %  'C:\Matlab\R2009b';    % 'D:\MATLAB\R2011b'; %
sBaseRunID     = 'LAST';    % ''; % 'LAST'; % '20100205_170006_125'; %
[pTask pInitWave] = setTaskParams('Schrodinger', '1D', sExampleName, sDataDirectory, sBaseRunID);
%
plotErrorNorm    = true;  % false; % 
plotSolutionNorm = false; % true;  % 
useRichardson    = false; % true;  % 
nRichardson      = 1;
nErrorFrequency  = 10;
sT   = 'FEM';   % 'FFDS';  % 
nN   = 2;       % 1/6;    %
nM_F = 1;       % 2.^(0:5);% 2.^(-4:0);   % 2.^0;     % 2.^(-2:1)
nJ_F = [1:0.5:6 6:3:60]; % [1:10 10:5:50 50:10:150];       % 2.^0;        % [1:2 4:5];   % 2.^(-5:0) % 2.^(-4:0)    %
switch sExampleName
    case 'EX01'
        nM = 3000;  % 12*16*16;  %
        nJ = 10;    % 100;   % 
        %
        sT_E = 'FEM';
        nN_E = 9;
        nM_E = 12*nM; % 10*nM;
        nJ_E = 60;
    case 'EX02m'
        nM = 1000;
        nJ = (7/12)*1200;
        %
        sT_E = 'FEM';
        nN_E = 9;
        nM_E = 10*nM;
        nJ_E = 60;
    case 'EX20'
        nM = 16000; % 500;       % 10*1600; % 
        nJ = 36;    % 36*5;      % 36*160;  %
        %
        sT_E = 'FEM';   % 'FFDS'
        nN_E = 5;       % 9;       % 1/12
        nM_E = 12*10*1600; % 12*nM;   % 10*nM;
        nJ_E = 5*36;    % 320*36
    case 'EX22'
        nM = 4000;   % 1000; %
        nJ = 30;     % 1200; %
        %
        sT_E = 'FEM';
        nN_E = 9;
        nM_E = 4*nM;
        nJ_E = 15*nJ;
    case 'EX24' % Moxley
        nM = 2*50000;  % Moxley: 50000;
        nJ = 10*16;    % Moxley: 100*16;
        %
        sT_E = 'FEM';
        nN_E = 9;
        nM_E = 4*50000;
        nJ_E = 10*16;
    case 'EX25' % Sullivan
        nM = 5*(2/5)*4*1000;
        nJ = 2*40;
        %
        sT_E = 'FEM';
        nN_E = 9;
        nM_E = 4*nM;
        nJ_E = 4*nJ;
end
P = []; s = 1;
for m=nM_F
    for j=nJ_F
        P(s, :) = [ j*nJ m*nM ];
        s = s + 1;
    end
end
if strcmp(sT, 'FFDS'), sMatrixMethod = 'TRISYS'; else sMatrixMethod = 'QR'; end
for pp=1:size(P, 1)
    fprintf('pp=%d (%.2f)\n', pp, pp/size(P, 1)*100);

    %%
    P(pp, 4:7 ) = [ max(rError.CAbs) max(rError.L2Abs) max(rError.CRel) max(rError.L2Rel) ];
    P(pp, 8:11) = [ rError.CAbs(end) rError.L2Abs(end) rError.CRel(end) rError.L2Rel(end) ];
    P(pp,   12) = nTimeSolution;
    disp(nTimeSolution/60)
end
csvwrite([ sFileName '.csv' ], P);