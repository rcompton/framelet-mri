%
% CS-framelet recon with general operator A
%
close all;clear all;

addpath('./Framelet/'); %Jai Feng Cai's Framelet library
addpath('./BregmanCookbook/'); %Jerome Gilles' Bregman library (has some framelet stuff)
addpath('./nufft_files'); %Fessler+Greengard+Lustig's nufft

%use these all the time...
vec = inline('reshape(x,[numel(x) 1])','x');
unvec = inline('reshape(x,[m n])','x','m','n');

stream = RandStream('mt19937ar');

%img = double(rgb2gray(imread('bouchard_mri_clean.png')));
%img = double(imread('phantom.gif'));
%load mri;
%img = double(D(:,:,1,19));
%img = double(rgb2gray(imread('cleanbrain.png')));
%img = double(imread('brainweb_t1.jpg'));

img = phantom('Modified Shepp-Logan',128);

%load raw' data.mat';
%img = abs(fftshift(fft2(brain)));

%resize to a nice square
n = 128;
img = imresize(img,[n n]);


%noise
%sig = .01;
%img = img + sig*randn(size(img));

%number of sample for compressed sensing

%for dsamp = 1:8
%    for ldex = 1:4

dsamp = 5;
num_samples = round(n*n/dsamp);
R = zeros(n,n);
R(sort(randsample(1:n^2,num_samples))) = 1;

%% Load in the compressed E

load onetwentyeightE.mat

B = imresize(B,[n*n, size(B,2)]);
C = imresize(C,[size(C,1), n*n]);

A = @(x) samplefun(R,B,C,x,0);
At = @(x) samplefun(R,B,C,x,1);

%check adjoint
xx = randn(n*n,1);
yy = randn(n*n,1);
fprintf('relative transpose error: %d\n',(dot(A(xx),yy) - dot(xx,At(yy)))/dot(A(xx),yy));

%% Now that we have a system operator, sample the data in k space

%create sample data
sig = .01;
f = A(vec(img + sig*randn(size(img))));

%not sure why Tom did this normalization. It works really well though.
normFactor = n/norm(f(:));
f = f*normFactor;
img = (img/norm(img))*norm(f);


%We minimize nu*|nabla u| + exci*|Du| st Au = f
nu = 1;
exci = .0;

%nu = 1*(ldex==1 || ldex==3 ||ldex==4);
%exci = 1*(ldex==2 || ldex==3) + .01*(ldex==4);

%the splitting parameters. There's probably another optimization to figure
%how to pick the best ones...
mu = .5;
lambda = .5;
gamma = 5;

%the number of times I'm willing to apply A
maxiters = 3000;

%watch a video as you do all this?
video = 1;
imagesc(img + sig*randn(size(img)));
colormap bone
figure()

%% Run the reconstruction
tic
[u, errors, res_errors, errors_per_breg, res_error_per_breg, iters ] = bregman_cs_framelet_2dv2(f, n, n, A, At, nu, exci, mu, lambda, gamma, maxiters, img, video);
timenoq = toc
%%

h = figure
imagesc(real(u));
colormap bone;
axis off;
%print(h,'-dpng',sprintf('recon_nu%i_exci%f_dsamp%i.png',[nu exci dsamp]));
%save(sprintf('recon_nu%i_exci%f_dsamp%i.png',[nu exci dsamp]))
    

