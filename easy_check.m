%?
clear all;

n=4^2;

G = dftmtx(n);
E = randn(n,3);
E = E*E';

A = E.*G;

x = rand(n,1);

imagesc(E'*E*x);
figure();
imagesc(abs(A'*A*x));