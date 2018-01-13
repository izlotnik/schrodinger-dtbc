function [ E_C, E_L, E_Cr, E_Lr ] = loadMaxError(sFileName)
%
load(sFileName);
%
s = [ iJ iN ] - 1; % s = [ size(e_1, 2) size(e_1, 3) ];
%
E_C = transpose(reshape(e_1, s)); E_Cr = transpose(reshape(e_1r, s));
E_L = transpose(reshape(e_2, s)); E_Lr = transpose(reshape(e_2r, s));
%
E_C(E_C>1) = NaN; E_Cr(E_Cr>1) = NaN; 
E_L(E_L>1) = NaN; E_Lr(E_Lr>1) = NaN;
%
E_C = [ transpose([ NaN nNc ]) [ nJc; E_C ] ]; E_Cr = [ transpose([ NaN nNc ]) [ nJc; E_Cr ] ];
E_L = [ transpose([ NaN nNc ]) [ nJc; E_L ] ]; E_Lr = [ transpose([ NaN nNc ]) [ nJc; E_Lr ] ];