function c = calcDot(a,b)
%calcDot  Vector dot product.
%   C = calcDot(A,B) returns the scalar product of the vectors A and B.
%   A and B must be vectors of the same length.  When A and B are both
%   column vectors, calcDot(A,B) is the same as A.*B.
%
% function is based on MATLAB's dot function
a = a(:);
c = sum(a.*b);