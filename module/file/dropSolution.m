clearContext; format short; echo on;
%%
sExampleName   = 'EX22a';               % 'EX01'; % 'EX23'; % 
sDataDirectory = 'C:\MATLAB\R2011b';    % 'C:\Matlab\R2009b';    % 'D:\MATLAB\R2011b'; % 
sBaseRunID     = '';                % ''; % 'LAST'; % '20100205_170006_125'; %
[ pTask pInitWave ] = setTaskParams('Schrodinger', '2D', sExampleName, sDataDirectory, sBaseRunID); pTaskE = pTask;
%%
P = [ ...
%    1*1200/2 1*(3/5)*1000/2 1*128/4
%        1*1200/4/10 1*(3/5)*1000/2 1*128/4
%         16*1200 16*1000
%        1*1200/2 1*(3/5)*1000/2 1*128/4
        1*1200/1 1*(3/5)*1000/1 1*64/1
%    4*1200/1 2*(3/5)*1000/1 1*128/1
%        1*1200/1 2*(3/5)*1000/1 1*128/1       
%        1*1200/1 4*(3/5)*1000/1 1*128/1       
    ];
% s = size(P, 1)+1;
% for j=2.^(-2:1)
%     for m=2.^(-2:1)
%         for k=2.^(-2:1)
%             P(s, :) = [ j*1200 m*(3/5)*1000 k*128 ];
%             s = s + 1;
%         end
%     end
% end
for pp=1:size(P, 1)
    pTask = pTaskE; pCalc = []; pFile = []; pVisual = []; pp
    [pCalc   pFile] = setCalcParams(pTask, {'FFDS'}, [0 P(pp, 1) P(pp, 2) NaN 0 0 P(pp, 3)], {'TRISYS','SSP','CONST'}); % 
    % [pCalc   pFile] = setCalcParams(pTask, {'FFDS'}, [0 P(pp, 1) P(pp, 2) NaN 1], {'TRISYS','SSP'}); % ,'CONST'
    rmdir(pCalc.sBaseDirectory, 's')
end