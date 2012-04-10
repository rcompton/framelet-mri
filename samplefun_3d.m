function [ y ] = samplefun_3d(R,B,C,x,t)
%sum Bl*G*Cl x
%C is a cell array of mxnxk arrays. The arrays correspond to the diag(C) that you
%see in the papers. ie each 3D mxnxk block is a C_l
%
%
%input x is a column vector (you need this for pcg etc.)
%
%output is a column vector
%t is for transpose

[m,n,k] = size(R);

x = reshape(x,[m n k]);

L = numel(B);

%scale = sqrt(m*n*k);
scale = sqrt(numel(find(R)));%??

y = zeros(m,n,k);
if t
    for l=1:L
        Bx = B{l}.*x;
        y = y + C{l}.*ifftn( R.*Bx )*scale;
    end
    
else
    for l=1:L
        Cx = C{l}.*x;
        y = y + B{l}.*R.*fftn( Cx )/scale;
    end
end

y = reshape(y,[numel(y) 1]);
end
