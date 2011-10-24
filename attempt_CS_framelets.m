%
% Test out if I can CS with framelets
%
close all;clear all;clc;

addpath('./Framelet/');
addpath('./BregmanCookbook/');
addpath(genpath('~./Dropbox/irt'));


vec = inline('reshape(x,[numel(x) 1])','x');
unvec = inline('reshape(x,[m n])','x','m','n');

%img = double(rgb2gray(imread('bouchard_mri_clean.png')));
img = double(imread('phantom.gif'));
%load mri;
%img = double(D(:,:,1,21));
%img = double(rgb2gray(imread('cleanbrain.png')));
%img = double(imread('brainweb_t1.jpg'));



%crop because I need it to be a square.
%n = min(size(img));
n = 64;
img = imresize(img,[n n]);

[m n] = size(img);
assert(m==n);

num_samples = round(m*n/2.5);

R = zeros(m,n);
R(randsample(1:m*n, num_samples)) = 1.0;
scale = sqrt(m*n); %only needed when using matlab FFT

%A = @(x) vec(R.*fft2(unvec(x,m,n))/scale);
%At = @(x) vec(ifft2(R.*unvec(x,m,n))*scale);

%These work on column vectors
 %L = 4;
 %B = ones(n,n*L)./L;
 %C = ones(n,n*L)./L;
 %A = @(x) samplefun(R,B,C,x,0);
 %At = @(x) samplefun(R,B,C,x,1);

 
 %using nufft
 %define the freq data locations
[k1pts k2pts] = meshgrid(1:m, 1:n);
k1pts = reshape(k1pts,numel(k1pts),1)*2*pi/m;
k2pts = reshape(k2pts,numel(k2pts),1)*2*pi/n;
omega = [k1pts k2pts];
omega = omega(sort(randsample(1:max(size(omega)),num_samples)),:);
 
%the number of neighbors to use when interpolating
j1 = 5;
j2 = 5;
%the fft sizes?
k1 = 2*m;
k2 = 2*n;

%make the nufft object
args = {omega, [m n], [j1 j2], [k1 k2]};
st = nufft_init(omega, [m n], [j1 j2], [k1 k2]);

scale_factor = sqrt(num_samples);

%input vectors and output vectors as needed by pcg. Note that the output of
%nufft_adj has m*n elements but nufft outputs num_samples x 1 column vecotr
%A = @(x) nufft(unvec(x,m,n),st)./scale_factor;
%At = @(x) vec(nufft_adj(x,st)./scale_factor);
L = 1;
B = cell(1,L);
for i=1:L
    B{i} = ones(num_samples,1)/L;
end
C = cell(1,L);
for i=1:L
    C{i} = ones(m*n,1)/L;
end
A = @(x) samplefun_nufft(st,B,C,x,m,n,0);
At = @(x) samplefun_nufft(st,B,C,x,m,n,1);

f = A(vec(img));
%not sure why Tom did this
normFactor = 1/norm(f(:)/size(R==1,1));
f = f*normFactor;

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
mu = 1;
lambda = 1;
gamma = 5;

if nu == 0
    lambda = 0;
end
if exci == 0
    gamma = 0;
end

%create the AtA operator. for some reason lk convolution works best.
lk = zeros(m,n);
lk(1,1) = 4;lk(1,2)=-1;lk(2,1)=-1;lk(m,1)=-1;lk(1,n)=-1;

%AtA = @(x) vec( mu*At(A( unvec(x,m,n) )) + lambda*ifft2(fft2(reshape(x,m,n)).*fft2(lk)) + gamma.*unvec(x,m,n) );

%using the samplefun
%AtA = @(x) mu*At(A( x )) + vec(lambda*ifft2(fft2(unvec(x,m,n)).*fft2(lk))) + gamma.*x ;
%samplev = @(x)samplefun(R,B,C,x,false);
AtA = @(x) mu*At(A(x)) + vec( lambda*ifft2(fft2(unvec(x,m,n)).*fft2(lk))) + gamma*x;

%set initial guesses
U = FraDecMultiLevel(zeros(m,n),D,n_level);
dw = U;
bw = U;

u = unvec(At(vec(f)),m,n);
%u = randn(m,n);

dx = zeros(size(u));
dy = zeros(size(u));
bx = zeros(size(u));
by = zeros(size(u));


%constrained, replace f with fl and update after each unconstrained
fl = f;
errorsr = [0];
errors = [0];
subplot(2,2,1);
imagesc(img);colormap hot;
subplot(2,2,2);
for l = 1:1000
    %unconstrained
    for k=1:3
        %update u
        rhsD = FraRecMultiLevel(SubFrameletArray(dw,bw),Dt,n_level);
        rhs = mu.*unvec(At(vec(fl)),m,n) + lambda.*Dxt(dx - bx) + lambda.*Dyt(dy - by) + gamma.*rhsD;
        
        
        [u,flag,reles,iter] = pcg(AtA,vec(rhs),1e-3,20);
        
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
    %fl = fl + f - unvec(A(vec(u)),m,n);
    fl = fl + f - A(vec(u));
    
    errorsr = [errorsr norm(A(vec(u)) - f,'fro')/norm(f,'fro')];
    errors = [errors norm(u/max(u(:)) - img/max(img(:)),'fro')];


    if randi(1)==1
        subplot(2,2,2)
        imagesc(real(u))
        
        subplot(6,2,7);
        plot(errors(end-length(errors)/2 : end));
        subplot(6,2,9);
        plot(errorsr(end-length(errorsr)/2 : end-1),'r.-');

        
        subplot(2,2,4);
        errorfig = abs(u/norm(u) - img/norm(img));
        imagesc(errorfig);
        
        colormap hot;
        pause(0.01);
    
        fprintf('error (Au -f): %f \n', errors(end));
    end
    
end

figure()
plot(errors)
