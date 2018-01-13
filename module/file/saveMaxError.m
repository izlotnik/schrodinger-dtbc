clearContext; format short; echo on;

%% J=60_M=300
% sMnemo = {'J=60_M=300', '20130407_203122_1D_EX01_MaxError_N=1,2,3,4,5,6,7,8,9,10_M=300_FREQ=1_C'};
% [ E_C, E_L, E_Cr, E_Lr ] = loadMaxError(['data\plots\EX01_' sMnemo{1} '_RICH(2)\' sMnemo{2}]);
% % sList = 'RICH(4)'; sMnemo{2} = '20130407_210337_1D_EX01_MaxError_N=1,2,3,4,5,6,7,8,9,10_M=300_FREQ=1_C';
% sList = 'RICH(3)'; sMnemo{2} = '20130407_211002_1D_EX01_MaxError_N=1,2,3,4,5,6,7,8,9,10_M=300_FREQ=1_C';
% [ E_C_R, E_L_R, E_Cr_R, E_Lr_R ] = loadMaxError(['data\plots\EX01_' sMnemo{1} '_' sList '\' sMnemo{2}]);

%% J=60_M=3000
sMnemo = {'J=60_M=3000', '20130406_203720_1D_EX01_MaxError_N=1,2,3,4,5,6,7,8,9,10_M=3000_FREQ=1_C'};
[ E_C, E_L, E_Cr, E_Lr ] = loadMaxError(['data\plots\EX01_' sMnemo{1} '_RICH(1)\' sMnemo{2}]);
% sList = 'RICH(4)'; sMnemo{2} = '20130406_203456_1D_EX01_MaxError_N=1,2,3,4,5,6,7,8,9,10_M=3000_FREQ=1_C';
% sList = 'RICH(3)'; sMnemo{2} = '20130407_212524_1D_EX01_MaxError_N=1,2,3,4,5,6,7,8,9,10_M=3000_FREQ=1_C';
sList = 'RICH(2)'; sMnemo{2} = '20130407_204625_1D_EX01_MaxError_N=1,2,3,4,5,6,7,8,9,10_M=3000_FREQ=1_C';
[ E_C_R, E_L_R, E_Cr_R, E_Lr_R ] = loadMaxError(['data\plots\EX01_' sMnemo{1} '_' sList '\' sMnemo{2}]);

sMnemo(3:5)={'A','N','AA'}; % J=60
sFile = ['data\plots\_1D_EX01_MaxError_N=1,2,3,4,5,6,7,8,9,10_' sMnemo{1} '_FREQ=1.xls'];

%% J=780_M=300
% sMnemo = {'J=780_M=300', '20130407_190428_1D_EX01_MaxError_N=1,2,3,4,5_M=300_FREQ=1_C'};
% [ E_C, E_L, E_Cr, E_Lr ] = loadMaxError(['data\plots\EX01_' sMnemo{1} '_RICH(2)\' sMnemo{2}]);
% % sList = 'RICH(4)'; sMnemo{2} = '20130407_183250_1D_EX01_MaxError_N=1,2,3,4,5_M=300_FREQ=1_C';
% sList = 'RICH(3)'; sMnemo{2} = '20130407_185842_1D_EX01_MaxError_N=1,2,3,4,5_M=300_FREQ=1_C';
% [ E_C_R, E_L_R, E_Cr_R, E_Lr_R ] = loadMaxError(['data\plots\EX01_' sMnemo{1} '_' sList '\' sMnemo{2}]);

%% J=780_M=3000
% sMnemo = {'J=780_M=3000', '20130407_172902_1D_EX01_MaxError_N=1,2,3,4,5_M=3000_FREQ=1_C'};
% [ E_C, E_L, E_Cr, E_Lr ] = loadMaxError(['data\plots\EX01_' sMnemo{1} '_RICH(1)\' sMnemo{2}]);
% % sList = 'RICH(4)'; sMnemo{2} = '20130407_142356_1D_EX01_MaxError_N=1,2,3,4,5_M=3000_FREQ=1_C';
% % sList = 'RICH(3)'; sMnemo{2} = '20130407_172344_1D_EX01_MaxError_N=1,2,3,4,5_M=3000_FREQ=1_C';
% sList = 'RICH(2)'; sMnemo{2} = '20130407_144546_1D_EX01_MaxError_N=1,2,3,4,5_M=3000_FREQ=1_C';
% [ E_C_R, E_L_R, E_Cr_R, E_Lr_R ] = loadMaxError(['data\plots\EX01_' sMnemo{1} '_' sList '\' sMnemo{2}]);
%
% sMnemo(3:5)={'A','AB','BC'}; % J=780
% sFile = ['data\plots\_1D_EX01_MaxError_N=1,2,3,4,5_' sMnemo{1} '_FREQ=1.xls'];
%
sIndex = num2str(0*(size(E_C, 1)+1)+1);
xlswrite(sFile, [ E_C(:, 1) [ E_C(1, 2:end); E_C(2:end, 2:end)./E_C_R(2:end, 2:end) ] ], sList, [sMnemo{3} sIndex]); xlswrite(sFile, {'C-abs'}, sList, [sMnemo{3} sIndex]);
xlswrite(sFile, E_C_R, sList, [sMnemo{4} sIndex]); xlswrite(sFile, {sList}, sList, [sMnemo{4} sIndex]);
xlswrite(sFile, E_C  , sList, [sMnemo{5} sIndex]); xlswrite(sFile, {'FEM'}, sList, [sMnemo{5} sIndex]);
%
sIndex = num2str(1*(size(E_C, 1)+1)+1);
xlswrite(sFile, [ E_L(:, 1) [ E_L(1, 2:end); E_L(2:end, 2:end)./E_L_R(2:end, 2:end) ] ], sList, [sMnemo{3} sIndex]); xlswrite(sFile, {'L2-abs'}, sList, [sMnemo{3} sIndex]);
xlswrite(sFile, E_L_R, sList, [sMnemo{4} sIndex]); xlswrite(sFile, {sList}, sList, [sMnemo{4} sIndex]);
xlswrite(sFile, E_L  , sList, [sMnemo{5} sIndex]); xlswrite(sFile, {'FEM'}, sList, [sMnemo{5} sIndex]);
%
sIndex = num2str(2*(size(E_C, 1)+1)+1);
xlswrite(sFile, [ E_Cr(:, 1) [ E_Cr(1, 2:end); E_Cr(2:end, 2:end)./E_Cr_R(2:end, 2:end) ] ], sList, [sMnemo{3} sIndex]); xlswrite(sFile, {'C-rel'}, sList, [sMnemo{3} sIndex]);
xlswrite(sFile, E_Cr_R, sList, [sMnemo{4} sIndex]); xlswrite(sFile, {sList}, sList, [sMnemo{4} sIndex]);
xlswrite(sFile, E_Cr  , sList, [sMnemo{5} sIndex]); xlswrite(sFile, {'FEM'}, sList, [sMnemo{5} sIndex]);
%
sIndex = num2str(3*(size(E_C, 1)+1)+1);
xlswrite(sFile, [ E_Lr(:, 1) [ E_Lr(1, 2:end); E_Lr(2:end, 2:end)./E_Lr_R(2:end, 2:end) ] ], sList, [sMnemo{3} sIndex]); xlswrite(sFile, {'L2-rel'}, sList, [sMnemo{3} sIndex]);
xlswrite(sFile, E_Lr_R, sList, [sMnemo{4} sIndex]); xlswrite(sFile, {sList}, sList, [sMnemo{4} sIndex]);
xlswrite(sFile, E_Lr  , sList, [sMnemo{5} sIndex]); xlswrite(sFile, {'FEM'}, sList, [sMnemo{5} sIndex]);