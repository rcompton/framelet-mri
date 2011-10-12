%
% Test out if I can CS with framelets
%
close all;clear all;clc;

addpath('./Framelet/');
addpath('./BregmanCookbook/');

vec = inline('reshape(x,[numel(x) 1])','x');
unvec = inline('reshape(x,[m n])','x','m','n');

%img = double(rgb2gray(imread('bouchard_mri_clean.png')));
img = double(imread('phantom.gif'));
%load mri;
%img = double(D(:,:,1,13));

[m n] = size(img);

num_samples = round(m*n/4.1);

R = zeros(m,n);
R(randsample(1:m*n, num_samples)) = 1.0;
scale = sqrt(m*n); %only needed when using matlab FFT

A = @(x) R.*fft2(x)/scale;
At = @(x) ifft2(R.*x)*scale;


f = A(img);

%1 is for Haar, 2 is for framelet, 3 is for cubic framelet
%Haar sucks. 2 works best, 3 is slower and worse.
[D,Dt]=GenerateFrameletFilter(2);

%more levels is better? I don't know. Experimentally, more levels is slow
%and worse quality. Hmm.
n_level = 1;


%We minimize nu*|nabla u| + exci*|Du| st Au = f
nu = 1;
exci = 1;

%mu is the exterior constraint split
%gamma is the framelet split
%lambda is the TV split
%I think they all be .5
mu = .5;
lambda = .5;
gamma = .5;

if nu == 0
    lambda = 0;
end
if exci == 0
    gamma = 0;
end

%L = laplacian(m*n);
lk = zeros(m,n);
lk(1,1) = 4;lk(1,2)=-1;lk(2,1)=-1;lk(m,1)=-1;lk(1,n)=-1;
AtA = @(x) vec( mu*At(A( unvec(x,m,n) )) + lambda*ifft2(fft2(reshape(x,m,n)).*fft2(lk)) + gamma.*unvec(x,m,n) );

%set initial guesses
U = FraDecMultiLevel(zeros(m,n),D,n_level);
dw = U;
bw = U;
u = At(f);
dx = zeros(size(u));
dy = zeros(size(u));
bx = zeros(size(u));
by = zeros(size(u));


%constrained, replace f with fl and update after each unconstrained
tol = 1e-4;
fl = f;
errors = [tol tol*3];
for l = 1:1000
    %unconstrained
    for k=1:randi(3)
        %update u
        rhsD = FraRecMultiLevel(SubFrameletArray(dw,bw),Dt,n_level);
        rhs = mu.*At(fl) + lambda.*Dxt(dx - bx) + lambda.*Dyt(dy - by) + gamma.*rhsD;
        
        [u,flag,reles,iter] = pcg(AtA,vec(rhs),1e-4,20);
        if randi(5)==randi(5)
            fprintf('                                 pcg error: %d, iter: %i \n', [reles iter]);
        end
        
        u = unvec(u,m,n);
        
        %update d's
        U = AddFrameletArray(FraDecMultiLevel(u,D,n_level),bw); %Du + b_k
        if gamma ~= 0 && exci ~= 0
            dw = ShrinkFramelet(U,1/(gamma/exci));
        end
        
        if lambda ~=0 && nu ~= 0
            [dx,dy] = shrink2( Dx(u)+bx, Dy(u)+by,1/(lambda/nu));
        end
        
        %update b's
        bw = SubFrameletArray(U,dw); %U-d_k
        bx = bx + (Dx(u) - dx);
        by = by + (Dy(u) - dy);
        
    end
    fl = fl + f - A(u);
    
    %surf(abs(u),'EdgeColor','None');
    %view(160,100)
    imagesc(real(u))
    colormap bone;
    pause(0.01);
    
    errors = [errors norm(u - img)/norm(img)];
    
    if randi(3)==randi(3)
        fprintf('error ( not residual): %f \n', errors(end));
    end
end

figure()
plot(errors)
