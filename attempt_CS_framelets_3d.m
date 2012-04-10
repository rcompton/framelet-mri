%
% CS-framelet recon with general operator A
%
close all;clear all;clc;

addpath('./Framelet/'); %Jai Feng Cai's Framelet library
addpath('./BregmanCookbook/'); %Jerome Gilles' Bregman library (has some framelet stuff)
addpath('./nufft_files'); %Fessler+Greengard+Lustig's nufft

%use these all the time...
vec = inline('reshape(x,[numel(x) 1])','x');
unvec = inline('reshape(x,[m, n, k])','x','m','n','k');

stream = RandStream('mt19937ar');

%img = double(rgb2gray(imread('bouchard_mri_clean.png')));
%img = double(imread('phantom.gif'));
%load mri;
%img = double(D(:,:,1,19));
%img = double(rgb2gray(imread('cleanbrain.png')));
%img = double(imread('brainweb_t1.jpg'));

% %resize to a nice square for a 2D experiment
% n = 64;
% m = 64;
% k = 1;
% img = imresize(img,[m n]);

%easy 3D
 img = zeros(32,32,32);
 img(12:20,8:15:9:20) = 1.0;
 img(20:30,8:15,3:7) = 2.0;
 img(15:24,20:30,18:27) = 3.8;
 img = img + shiftdim(img,2);


load littlebabymri.mat;
img = dee;

[m n k] = size(img);

%number of sample for compressed sense
num_samples = round(m*n*k/1.5);


%% Create sample operator A
% % A is a function handle, you'll need to make At as well.
% % Note that input vectors and output vectors as needed by pcg. Note that the output of
% % nufft_adj has m*n elements but nufft outputs num_samples x 1 column vecotr
% 
% %Create B and C, this is a well understood research problem and I'm not going to
% %bother with making a better interpolator
% L = 3;
% T = .2;
% mt = num_samples;
% t = linspace(0.1,T,mt); %number of time points in the time discreization of continuous time
% tl = linspace(0.1,T,L); %number of time points in to approximate at
% 
% %load the other guy's field map
% load fmap_test.mat
% fmap = imresize(fmap,size(img));
% fmap = rot90(fmap,-1);
% 
% w = fmap/norm(fmap(:));
% 
% C = cell(1,L);
% for l=1:L
%     C{l} = vec(exp(1j*w*tl(l)));
% end
% 
% B = cell(1,L);
% for l=1:L
%     B{l} = zeros(num_samples,1);
% end
% for l=1:L-1
%     for i=1:num_samples
%         if (tl(l) <= t(i)) && (t(i) <= tl(l+1))
%             %B(i,(l-1)*N+1:l*N) = 1 - (t(i) - tl(l))/(tl(l+1) - tl(l));
%             %B(i,l+1) = (t(i) - tl(l))/(tl(l+1) - tl(l));
%             
%             B{l}(i) = 1 - (t(i) - tl(l))/(tl(l+1) - tl(l));
%             B{l+1}(i) = (t(i) - tl(l))/(tl(l+1) - tl(l));
%         end
%     end
% end
% 
% %create sample operators
% A = @(x) samplefun_nufft(st,B,C,x,m,n,0);
% At = @(x) samplefun_nufft(st,B,C,x,m,n,1);
% 
% % 
% %without nufft library
% %  fprintf('NO NUFFT!\n');
R = zeros(m,n,k);
sample_inds = randsample(stream, 1:m*n*k, num_samples);
R(sample_inds) = 1.0;
% % 
%  B2 = zeros([size(R,1) L*size(R,2)]);
%  C2 = zeros([size(R,1) L*size(R,2)])/L;
%  for l=1:L
%      Bee = zeros(size(R));
%      Bee(sample_inds) = B{l};
%      B2(:, (l-1)*n +1 : l*n) = Bee;
%      
%      C2(:, (l-1)*n +1 : l*n) = reshape(C{l},[m n]);
% 
%  end
%      
%  A = @(x) samplefun(R,B2,C2,x,0);
%  At = @(x) samplefun(R,B2,C2,x,1);

fprintf('NO CORRECTION!\n');

%compresses sensing operator scale factor
scale_factor = sqrt(num_samples);

A = @(x) reshape( R.*fftn( reshape(x,[m n k]) )/scale_factor, [numel(x) 1]);
At = @(x) reshape( ifftn( R.*reshape(x,[m n k]) )*scale_factor, [numel(x) 1]);


%% Now that we have a system operator, sample the data in k space

%create sample data
%f = single(A(vec(img)));
f = A(reshape(img,[numel(img) 1]));

%not sure why Tom did this normalization. It works really well though.
normFactor = 1/norm(f(:)/m);
f = f*normFactor;


%We minimize nu*|nabla u| + exci*|Du| st Au = f
nu = 1;
exci = 1;

%the splitting parameters. There's probably another optimization to figure
%how to pick the best ones...
mu = 1;
lambda = 0.5;
gammah = 5;
%gammah = 0;

%the number of times I'm willing to apply A
maxiters = 33231;

%watch a video as you do all this?
video = 1;

%% Run the reconstruction
tic
[u, errors, res_errors, errors_per_breg, res_error_per_breg, iters ] = bregman_cs_framelet_3d(f, m, n, k, A, At, nu, exci, mu, lambda, gammah, maxiters, img, video);
toc

