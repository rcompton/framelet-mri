function Rec=FraRecMultiLevel(C,R,L)

% function Rec=FraRecMultiLevel(C,R,L)
% This function implement framelet reconstruction up to level L.
% C ==== the data to be reconstructed, which are in cells in C{i,j} with
% C{1,1} being a cell.
% R ==== is the reconstruction filter in 1D. In 2D, it is generated by tensor
% product. The filter D must be symmetric or anti-symmetric, which
% are indicated by 's' and 'a' respectively in the last cell of R.
% L ==== is the level of the decomposition.
% Rec ==== is the reconstructed data.

% Written by Jian-Feng Cai.
% email: tslcaij@nus.edu.sg
%
% See FraDec FraRec FraDecMultiLevel GenerateFrameletFilter

nR=length(R);

for k=L:-1:2
    C{k-1}{1,1}=FraRec(C{k},R,k);
end
Rec=FraRec(C{1},R,1);