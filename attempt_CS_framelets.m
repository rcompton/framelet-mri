%
% Test out if I can CS with framelets
%
close all;clear all;clc;

addpath('./Framelet/');
addpath('./BregmanCookbook/');

vec = inline('reshape(x,[numel(x) 1])','x');
unvec = inline('reshape(x,[m n])','x','m','n');

%img = double(rgb2gray(imread('../bouchard_mri_clean.png')));
img = double(imread('phantom.gif'));
[m n] = size(img);

num_samples = round(m*n/4);

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
[D,Dt]=GenerateFrameletFilter(2);

%more levels is better? I don't know
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

n_outer = 1250;
n_inner = 1;
errors = zeros(n_outer, 1);

%constrained, replace f with fl and update after each unconstrained

fl = f;
for l=1:n_outer

    %unconstrained
    for k=1:randi(20),
        %update u
        rhs = FraRecMultiLevel(SubFrameletArray(dk,bk),Dt,n_level);
        rhs = mu*At(fl) + lambda.*rhs;
        
        [u,~] = pcg(AtA,vec(rhs),1e-3,20);
        u = unvec(u,m,n);
               
        %update d
        U = AddFrameletArray(FraDecMultiLevel(u,D,n_level),bk); %Du + b_k
        dk = ShrinkFramelet(U,1/lambda);
        
        %update b
        bk = SubFrameletArray(U,dk); %U-d_k
        
        
    end
    fl = fl + f - A(u);
    
    imagesc(real(u));
    colormap bone;
    pause(0.015);
    
    errors(l) = norm(u-img,'fro');
    fprintf('error (not residual): %f \n', errors(l));

end

figure()
plot(errors)
