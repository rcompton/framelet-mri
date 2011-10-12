%
% Test out if I can CS with framelets
%
close all;clear all;clc;

addpath('./Framelet/');
addpath('./BregmanCookbook/');

vec = inline('reshape(x,[numel(x) 1])','x');
unvec = inline('reshape(x,[m n])','x','m','n');

img = double(rgb2gray(imread('bouchard_mri_clean.png')));
%img = double(imread('phantom.gif'));
%load mri;
%img = double(D(:,:,1,13));

[m n] = size(img);

num_samples = round(m*n/pi);

R = zeros(m,n);
R(randsample(1:m*n, num_samples)) = 1.0;
A = @(x) R.*fft2(x);
At = @(x) ifft2(R.*x);
%
% Am = randn(num_samples, m*n);
% A = @(x) Am*vec(x);
% At = @(x) unvec( Am'*x, m, n);

f = A(img);

%1 is for Haar, 2 is for framelet, 3 is for cubic framelet
%Haar sucks. 2 works best, 3 is slower and worse.
[D,Dt]=GenerateFrameletFilter(2);

%more levels is better? I don't know. Experimentally, more levels is slow
%and worse quality. Hmm.
n_level = 1;

%how to framelet decompose and recompose. The decomposed coeff are in a
%cell or somthing
%Dimg = FraDecMultiLevel(img,D,n_level);
%DtDimg = FraRecMultiLevel(Dimg,Dt,n_level);



%AtA = @(x) vec( At(A( unvec(x,m,n) )) + lambda.*unvec(x,m,n) );

mu = .5;
lambda = .5;

% a = .3;
% b = .6;    
% mu = (b-a)*rand() + a
% lambda = (b-a)*rand() + a

AtA = @(x) vec( mu*At(A( unvec(x,m,n) )) + lambda.*unvec(x,m,n) );
%AtA = mu*(Am'*Am) + lambda*eye(m*n);

%
% %wtf?
U = FraDecMultiLevel(zeros(m,n),D,n_level);
dk = U;
bk = U;
u = norm(f).*randn(size(f));


errors = [1];

%constrained, replace f with fl and update after each unconstrained
tol = 1e-2;
fl = f;
while errors(end) > tol

    %unconstrained
    for k=1:randi(1)
        %update u
        rhs = FraRecMultiLevel(SubFrameletArray(dk,bk),Dt,n_level);
        rhs = mu*At(fl) + lambda.*rhs;
        
        [u,flag,reles,iter] = pcg(AtA,vec(rhs),1e-4,20);
        if randi(5)==randi(5)
            fprintf('                                 pcg error: %d, iter: %i \n', [reles iter]);
        end

        u = unvec(u,m,n);
               
        %update d
        U = AddFrameletArray(FraDecMultiLevel(u,D,n_level),bk); %Du + b_k
        dk = ShrinkFramelet(U,1/lambda);
        
        %update b
        bk = SubFrameletArray(U,dk); %U-d_k
        
        
    end
    fl = fl + f - A(u);
    
    %surf(abs(u),'EdgeColor','None');
    %view(160,100)
    imagesc(abs(u))
    colormap bone;
    pause(0.015);
    
    errors = [errors norm(u - img)/norm(img)];
    
    if randi(3)==randi(3)
        fprintf('error ( not residual): %f \n', errors(end));
    end
end

figure()
plot(errors)
