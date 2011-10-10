function u=Framelet_NB_Deconvolution2(f,A,mu,lambda,delta,Niter,frame,NLevel,typeKernel)

%==========================================================================
%
% u=Framelet_NB_Deconvolution2(f,A,mu,lambda,Niter,frame,NLevel,typeKernel)
%
% This function performs the Nonblind Framelet deconvolution by
% the synthesis approach.
% -v 1.2: 07/18/2011
% -v 1.0: 05/05/2011
%
% Parameters:
% f: input blurred image
% A: input blur kernel
% mu,lambda: regularization parameters
% delta: gradient descent speed
% Niter: number of iterations
% frame: type of used Framelet (0=Haar, 1=Piecewise Linear
%        Framelet, 2=Piecewise Cubic Framelet)
% NLevel: number of scale used in the Framelet decomposition
% typeKernel: 0 <=> A is considered in the spatial domain
%             1 <=> A is considered in the Fourier domain and must be of
%             the same size of the image with frequency (0,0) at location
%             (1,1)
%
% NB: mu=1000,lambda=10,delta=20,Niter=50 with Haar and NLevel=3
% perform well for a 3x3 disc kernel.
%
% Author: Jerome Gilles
% Institution: UCLA - Math Department
% email: jegilles@math.ucla.edu
%
%===============================================================

%Generation of Framelet filters
[D,R]=GenerateFrameletFilter(frame);
nD=length(D);

%structures initialization
[M,N]=size(f);
U=FraDecMultiLevel(zeros(M,N),D,NLevel);
d=U;
c=U;

%Fourier mask + initialization of Fourier constant quantities
if typeKernel==0
    Mask=zeros(M,N);
    [H,L]=size(A);
    Mask([end+1-floor(H/2):end,1:ceil(H/2)],[end+1-floor(L/2):end,1:ceil(L/2)]) = A;
    FMask=fft2(Mask);
else
    FMask = A;
end

P=(abs(FMask).^2+lambda).^-1;
FW=delta*conj(FMask).*P.*FMask;
FF=delta*conj(FMask).*P.*fft2(f);

K=0;

while K<Niter,
    K=K+1
    %update c
    tx=FW.*fft2(FraRecMultiLevel(d,R,NLevel))-FF;
    c=SubFrameletArray(c,FraDecMultiLevel(real(ifft2(tx)),D,NLevel));
    
    %update d
    d=ShrinkFramelet(c,1/mu);
end

tu=FraRecMultiLevel(d,R,NLevel);
u=bilateralFilter(tu);
