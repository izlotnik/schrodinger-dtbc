function b = calcIDST(a)
%% calcIDST  
%
% function is based on MATLAB's idst function
n = size(a, 1) + 1;
b = ( 2 / n ) * calcDST(a);