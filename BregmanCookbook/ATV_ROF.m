function u=ATV_ROF(f,mu,lambda,Niter)

%============================================
% function u=ATV_ROF(f,mu,lambda,eps)
%
% Anisotropic ROF denoising
% Version:
% -v1.0 - 04/04/2011
%
% This function performs the minimization of
% u=arg min |D_x u|+|D_y u|+0.5*mu*||u-f||_2^2
%
% by the Split Bregman Iteration
%
% f = noisy image
% mu = regularization parameter
% lambda = fidelity factor for the split variables
% Niter = number of iteration
%
% Author: Jerome Gilles
% Institution: UCLA - Math Department
% email: jegilles@math.ucla.edu
%
% Note: mu=10, lambda=5 and Niter=5 perform well
%
%============================================

[M,N]=size(f);

f=double(f);
dx=zeros(M,N);
dy=zeros(M,N);
bx=zeros(M,N);
by=zeros(M,N);
u=f;
Z=zeros(M,N);
up=ones(M,N);

a=lambda/(mu+8*lambda);
b=mu/lambda;
K=0;

while K<Niter,
    K=K+1;
    up=u;
    tx=dx-bx;
    ty=dy-by;
   
    %Update u
    u=a*(u(:,[1,1:N-1])+u(:,[2:N,N])+u([1,1:M-1],:)+u([2:M,M],:)+4*u+b*f-(tx-tx(:,[1,1:N-1])+ty-ty([1,1:M-1],:)));
    ux=u(:,[2:N,N])-u;
    uy=u([2:M,M],:)-u;
        
    %Update dx
    tmpx=ux+bx;
    dx=sign(tmpx).*max(Z,abs(tmpx)-1/lambda);
    
    %Update dy
    tmpy=uy+by;
    dy=sign(tmpy).*max(Z,abs(tmpy)-1/lambda);
   
    %Update bx and by 
    bx=tmpx-dx;
    by=tmpy-dy;
end
