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

%2D datasets
%img = double(rgb2gray(imread('bouchard_mri_clean.png')));
%img = double(imread('phantom.gif'));
%load mri;
%img = double(D(:,:,1,19));
%img = double(rgb2gray(imread('cleanbrain.png')));
%img = double(imread('brainweb_t1.jpg'));

%img = imresize(img,[128,128]);

%easy 3D
img = zeros(32,32,32);
img(12:20,8:15:9:20) = 1.0;
img(20:30,8:15,3:7) = 2.0;
img(15:24,20:30,18:27) = 3.8;
img = img + shiftdim(img,2);


%load littlebabymri.mat;
%img = dee;

load twentytwomri128.mat
img = deenoise;

%load babymri64.mat
%img = dee64;

%img = double(img);
[m n k] = size(img);

%number of sample for compressed sense
num_samples = round(m*n*k/2.5);


%% Create sample operator A
% % A is a function handle, you'll need to make At as well.
% % Note that input vectors and output vectors as needed by pcg. Note that the output of
% % nufft_adj has m*n elements but nufft outputs num_samples x 1 column vecotr
 
% %Create B and C, this is a well understood research problem and I'm not going to
% %bother with making a better interpolator
L = 3;
T = .2;
mt = num_samples;
t = linspace(0.1,T,mt); %number of time points in the time discreization of continuous time
tl = linspace(0.1,T,L); %number of time points in to approximate at

% %load the other guy's field map
 load fmap_test.mat
 fmap = imresize(fmap,[m n]);
 fmap = rot90(fmap,-1);
 
efmap = zeros(m,n,k);
for kay=1:k
    efmap(:,:,kay) = fmap*( -(kay/k - 1)*(kay/k) );
end
fmap = efmap + shiftdim(efmap,2);

w = (fmap-min(fmap(:)))/max(fmap(:));
 
C = cell(1,L);
for l=1:L
     C{l} = exp(1j*w*tl(l));
end
B = cell(1,L);
for l=1:L
     B{l} = zeros(num_samples,1);
end
for l=1:L-1
     for i=1:num_samples
         if (tl(l) <= t(i)) && (t(i) <= tl(l+1))
             B{l}(i) = 1 - (t(i) - tl(l))/(tl(l+1) - tl(l));
             B{l+1}(i) = (t(i) - tl(l))/(tl(l+1) - tl(l));
        end
    end
end

R = zeros(m,n,k);
sample_inds = randsample(stream, 1:m*n*k, num_samples);
R(sample_inds) = 1.0;

%adjust away from that nufft format..
for l=1:L
     Bee = zeros(size(R));
     Bee(sample_inds) = B{l};
     B{l} = Bee;     
end


%fprintf('NO CORRECTION!\n');
%scale_factor = sqrt(num_samples);
%scale_factor = sqrt(m*n*k);
%A = @(x) reshape( R.*fftn( reshape(x,[m n k]) )/scale_factor, [numel(x) 1]);
%At = @(x) reshape( ifftn( R.*reshape(x,[m n k]) )*scale_factor, [numel(x) 1]);

%fprintf('3D crazy sample function!\n');
A = @(x) samplefun_3d(R,B,C,x,0);
At = @(x) samplefun_3d(R,B,C,x,1);

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
lambda = 1;
gammah = 1;
%gammah = 0;

%the number of times I'm willing to apply A
maxiters = 2000;

%number of bregs
max_bregs = 70;

%watch a video as you do all this?
video = 0;

%% Run the reconstruction
tic
%[u, errors, res_errors, errors_per_breg, res_error_per_breg, iters ] = bregman_cs_framelet_3d(f, m, n, k, A, At, nu, exci, mu, lambda, gammah, maxiters, img, video);
%[u, errors, res_errors, errors_per_breg, res_error_per_breg, iters ] = bregman_cs_framelet(f, m, n, A, At, nu, exci, mu, lambda, gammah, maxiters, img, video);

%bregman_cs_framelet
%
%We minimize nu*|nabla u| + exci*|Fu| st Au = f
%F is a framelet transform
%
% A -- a function handle to the sampling operator
% At -- A^t
%
% nu -- tune the TV strength
% exci -- tune the framelet strength
%
% mu -- the exterior constraint split
% gamma -- the framelet split
% lambda -- the TV split
% ie mu,lambda, and gamma are parameters for
% mu*At(A(x)) + lambda*laplacian + gamma*x
% which we solve with pcg
%
% maxiter -- how many times we can apply A before stopping the algorithm
%
% img is the exact image, we only use it to update the errors as you go.
%

%don't know if this exists already or not...
norm3 = @(x) sqrt(trapz(trapz(trapz(abs(x).^2))));


%1 is for Haar, 2 is for framelet, 3 is for cubic framelet
%Haar sucks. 2 works best, 3 is slower and worse.
%[D,Dt]=GenerateFrameletFilter(2);

%more levels is better? I don't know. Experimentally, more levels is slow
%and worse quality. Hmmm.
%n_level = 1;

%check that we chose sane parameters
if nu == 0
    assert(lambda == 0);
end
if exci == 0
    assert(gammah == 0);
end

%create the AtA operator. for some reason lk convolution like this (not conv2) works best.
%lk = zeros(m,n);
%lk(1,1) = 4;lk(1,2)=-1;lk(2,1)=-1;lk(m,1)=-1;lk(1,n)=-1;
%AtA = @(x) mu*At(A(x)) + lambda*reshape(ifft2(fft2(reshape(x,[m n])).*fft2(lk)), [m*n 1]) + gammah*x;

%2D
%Lap = laplacian([m n], {'P','P'}); %negative laplacian
%AtA = @(x) mu*At(A(x)) + lambda*Lap*x + gammah*x;


%3D, still seems to work for 2D data
fprintf('3D Laplacian matrix with 7-point stencil\n');
Lap = laplacian([m n k], {'P','P','P'});
AtA = @(x) mu*At(A(x)) + lambda*Lap*x + gammah*x;

%without correction or pcg
%     uker = zeros(m,n,k);
%     uker(1,1,1) = 6;
%     uker(2,1,1) = -1; uker(1,2,1) = -1; 
%     if k>1
%        uker(1,1,2) = -1;
%        uker(1,1,k) = -1;
%     end
%     uker(m,1,1) = -1; uker(1,n,1)  =-1; 
%     uker = mu*(conj(R).*R)+lambda*fftn(uker)+gammah;

%initial guess.
u = reshape(At(f),[m n k]);

%initial splitting vectors
dx = zeros(size(u));
dy = zeros(size(u));
dz = zeros(size(u));
bx = zeros(size(u));
by = zeros(size(u));
bz = zeros(size(u));

%constrained, replace f with fl and update after each unconstrained
fl = f;
res_errors = [1];
errors = [];
errors_per_breg = [];
res_errors_per_breg = [];
iters = [];

alphaflag=1;

% start the optimization.
% the outer loop is constraint enforcement
ell = 0;
while(sum(iters) < maxiters && res_errors(end)>1e-5 && ell < max_bregs)
    ell = ell+1;
    %unconstrained step, 5 is overkill.
    for kay=1:3
        %update u
        rhs = mu.*reshape(At(reshape(fl, [numel(fl) 1])),[m n k]);
        rhs = rhs + lambda*(Deriv(dx - bx,1,true) + Deriv(dy - by,2,true) + Deriv(dz - bz,3,true));
        rhs = rhs + gammah.*u; %Haar wavelet trick like in Tom's paper
        
        %this is where all the computation happens
        macpcg = 150;
        [u,flagcg,~,iterc] = pcg(AtA, reshape(rhs,[numel(rhs) 1]), 1e-3, macpcg);
        
        %Tom's way
  %       iterc = 1;
  %       u = ifftn(fftn(rhs)./uker);

        if flagcg==0
            fprintf('pcg successful, iterc: %i\n',iterc);        
            iters = [iters iterc];

        else
            fprintf('pcg terrible: %i\n',macpcg);
            iters = [iters macpcg];
        end
        
        u = reshape(u,[m n k]);
        
        %update d's
        %U = AddFrameletArray(FraDecMultiLevel(u,D,n_level),bw); %Du + b_k
        %if gammah ~= 0 && exci ~= 0
        %    dw = ShrinkFramelet(U,1/(gamma/exci));
        %end
        
        if lambda ~=0 && nu ~= 0
            [dx,dy,dz] = shrink3( Deriv(u,1,false)+bx, Deriv(u,2,false)+by, Deriv(u,3,false)+bz, 1/(lambda/nu));
        end
       
        %update b's
        %bw = SubFrameletArray(U,dw); %U-d_k
        bx = bx + (Deriv(u,1,false) - dx);
        by = by + (Deriv(u,2,false) - dy);
        bz = bz + (Deriv(u,3,false) - dz);
        
        res_errors = [res_errors norm(A(reshape(u, [m*n*k 1])) - f,'fro')/norm(f,'fro')];
        errors = [errors norm3(u./max(u(:)) - img./max(img(:)))];
        
    end
    fl = fl + f - A(reshape(u, [m*n*k 1]));
    res_errors_per_breg = [res_errors_per_breg norm(A(reshape(u, [m*n*k 1])) - f,'fro')/norm(f,'fro')];
    errors_per_breg = [errors_per_breg norm3(u./max(u(:)) - img./max(img(:)))];
    
    if (video && randi(5) == 3)
        if(ndims(u)==3)
        %%
        hh = figure;
        set(hh,'Visible','off');
        
        whitebg('k');       
        h = vol3d('cdata',abs(u),'texture','3D');
        %view(3);
        view(magicview);
        axis tight;  daspect([1 1 .4])
        alphamap('rampup');
        %colormap bone;
        %%
        %alphamap(.05 .* alphamap);
            
        fnameh = sprintf('3drecon_step%i',ell);
        print(hh,'-dpng',fnameh);
        
        close(hh);
        %pause(.03);
        else if (ndims(u)==2)
                subplot(2,2,1);
                imagesc(real(img));
                title('exact image');
        
                subplot(2,2,2);
                imagesc(real(u));
                title('current reconstruction');
        
                subplot(2,2,4);
                errorfig = abs(u/max(u(:)) - img/max(img(:)));
                imagesc(errorfig);
                title('reconstruction error');
        
                subplot(6,2,7);
                semilogy(errors_per_breg);
                title('image domain error per breg');
                subplot(6,2,11);
                semilogy(res_errors_per_breg(2:end),'r.-');
                title('residual error per breg');
        
                colormap hot;
                pause(0.03);
        
            end
        end
    end
    fprintf('breg step ell = %i error (u-img): %f Ax iter numer: %i \n', [ell errors(end) sum(iters)]);
    fprintf('breg step ell = %i res error : %f Ax iter numer: %i \n---\n', [ell res_errors(end) sum(iters)]);
    
end

toc

