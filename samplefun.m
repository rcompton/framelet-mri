function [ y ] = samplefun(R,B,C,x,t)
%sum Bl*G*Cl x
%C is a M by NL matrix with blocks that correspond to the diag(C) that you
%see in the papers
%x is a column vector (you need this for pcg etc.)
%output is a column vector
%t is for transpose

[m nl] = size(B);
[mn one] = size(x);

n = mn/m;

x = reshape(x,m,n);

L = nl/n;

scale = sqrt(m*n);

y = zeros(m,n);
if t
    for l=1:L
        Bx = B(:, (l-1)*n +1 : l*n).*x;
        y = y + C(:,(l-1)*n +1 : l*n).*ifft2( R.*Bx )*scale;
    end
    
else
    for l=1:L
        Cx = C(:, (l-1)*n +1 : l*n).*x;
        y = y + B(:,(l-1)*n +1 : l*n).*R.*fft2( Cx )/scale;
    end
end

y = reshape(y,numel(y),1);
end

