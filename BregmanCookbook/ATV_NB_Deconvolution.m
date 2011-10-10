function u=ATV_NB_Deconvolution(f,A,mu,lambda,Niter,typeKernel)

%============================================
% function u=ATV_NB_Deconvolution(f,A,mu,lambda,Niter)
%
% Anisotropic TV Non-Blind Deconvolution
% Version:
% -v1.2 - 07/18/2011
% -v1.0 - 04/05/2011
%
% This function performs the minimization of
% u=arg min |D_x u|+|D_y u|+0.5*mu*||A°u-f||_2^2
%
% by Split Bregman Iteration
%
% f = noisy image
% A = convolution kernel
% mu = regularization parameter
% lambda = fidelity factor for the split variables
% Niter = number of iterations
% typeKernel: 0 <=> A is considered in the spatial domain
%             1 <=> A is considered in the Fourier domain and must be of
%             the same size of the image with frequency (0,0) at location
%             (1,1)
%
% Author: Jerome Gilles
% Institution: UCLA - Math Department
% email: jegilles@math.ucla.edu
%
% Note: for a gaussian blur with sigma=1,
%  mu=500, lambda=1 and Niter=5 perform well
%
%============================================

[M,N]=size(f);

%Structures and constants initialization
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

%Fourier mask + initialization of Fourier constant quantities
if typeKernel==0
    Mask=zeros(M,N);
    [H,L]=size(A);
    Mask([end+1-floor(H/2):end,1:ceil(H/2)],[end+1-floor(L/2):end,1:ceil(L/2)]) = A;
    FMask=fft2(Mask);
else
    FMask = A;
end

%Fourier Laplacian mask initialization
D = zeros(M,N);
D([end,1,2],[end,1,2]) = [0,1,0;1,-4,1;0,1,0];
FD=fft2(D);

%Fourier constant initialization
FW=((mu/lambda)*abs(FMask).^2-real(FD)).^-1;
FF=(mu/lambda)*conj(FMask).*fft2(f);

%Bregman iterations
while K<Niter,
    K=K+1;
    up=u;
    tx=dx-bx;
    ty=dy-by;

    %Update u
    u=real(ifft2(FW.*(FF-fft2(tx-tx(:,[1,1:N-1])+ty-ty([1,1:M-1],:)))));
    ux=u-u(:,[1,1:N-1]);
    uy=u-u([1,1:M-1],:);

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
