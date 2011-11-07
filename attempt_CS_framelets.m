%
% Test out if I can CS with framelets
%
close all;clear all;

addpath('./Framelet/');
addpath('./BregmanCookbook/');
addpath(genpath('~./Dropbox/irt'));

stream = RandStream('mt19937ar');


vec = inline('reshape(x,[numel(x) 1])','x');
unvec = inline('reshape(x,[m n])','x','m','n');

%img = double(rgb2gray(imread('bouchard_mri_clean.png')));
img = double(imread('phantom.gif'));
%load mri;
%img = double(D(:,:,1,21));
%img = double(rgb2gray(imread('cleanbrain.png')));
%img = double(imread('brainweb_t1.jpg'));

%resize to a nice square
n = 64;
img = imresize(img,[n n]);
m=n;

%number of sample for compressed sense
num_samples = round(m*n/3);

%the downsample operator matrix, could be done away with now that I have
%nufft library but it works and was easy
R = zeros(m,n);
R(randsample(stream, 1:m*n, num_samples)) = 1.0;

%normalizations are never consistent...
scale = sqrt(m*n);

 
%using nufft
%define the freq data locations
[k1pts k2pts] = meshgrid(1:m, 1:n);
k1pts = reshape(k1pts,numel(k1pts),1)*2*pi/m;
k2pts = reshape(k2pts,numel(k2pts),1)*2*pi/n;
omega = [k1pts k2pts];

%bias towards low frequency samples (sucks)
%samp_coords = zeros(1,m*n);
%j=1;
%i=m*n-1;
%while(j<=num_samples && i>1)
%    y = rand();
%    if y > i/(m*n)
%        samp_coords(i) = i;
%        j = j+1;
%    end
%    i = i-1;
%end
%samp_coords = find(samp_coords);

%uniformly sample
samp_coords = sort(randsample(stream, 1:max(size(omega)),num_samples));

omega = omega(samp_coords,:);

%%

%the number of neighbors to use when interpolating
j1 = 11;
j2 = 11;
%the fft sizes? 2*n works nice.
k1 = 2*m;
k2 = 2*n;

%make the nufft object
args = {omega, [m n], [j1 j2], [k1 k2]};
st = nufft_init(omega, [m n], [j1 j2], [k1 k2]);

%nufft is like fftw in that there is no normalization factor
scale_factor = sqrt(num_samples);

%Create sample operator
%input vectors and output vectors as needed by pcg. Note that the output of
%nufft_adj has m*n elements but nufft outputs num_samples x 1 column vecotr

% do nothing correction
% L = 6;
% B = cell(1,L);
% for i=1:L
%     B{i} = randn(num_samples,1)/sqrt(L);
% end
% C = cell(1,L);
% for i=1:L
%     C{i} = randn(m*n,1)/sqrt(L);
% end

%Create B and C, this is a solved research problem and I'm not going to
%bother with making a better interpolator
L = 3;
T = .2;
mt = num_samples;
t = linspace(0.1,T,mt); %number of time points in the time discreization of continuous time
tl = linspace(0.1,T,L); %number of time points in to approximate at

%load the other guys field map
load fmap_test.mat
fmap = imresize(fmap,size(img));

w = fmap/norm(fmap);
%w = peaks(n)/norm(peaks(n));
%w = img./norm(img(:));


C = cell(1,L);
for l=1:L
    C{l} = vec(exp(1j*w*tl(l)));
end

B = cell(1,L);
for l=1:L
    B{l} = zeros(num_samples,1);
end
for l=1:L-1
     for i=1:num_samples
         if (tl(l) <= t(i)) && (t(i) <= tl(l+1))
             %B(i,(l-1)*N+1:l*N) = 1 - (t(i) - tl(l))/(tl(l+1) - tl(l));
             %B(i,l+1) = (t(i) - tl(l))/(tl(l+1) - tl(l));
             
             B{l}(i) = 1 - (t(i) - tl(l))/(tl(l+1) - tl(l));
             B{l+1}(i) = (t(i) - tl(l))/(tl(l+1) - tl(l));
         end
     end
end

 %create sample operators
A = @(x) samplefun_nufft(st,B,C,x,m,n,0);
At = @(x) samplefun_nufft(st,B,C,x,m,n,1);

%create sample data
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
AtA = @(x) mu*At(A(x)) + lambda*vec(ifft2(fft2(unvec(x,m,n)).*fft2(lk))) + gamma*x;

%set initial guesses
U = FraDecMultiLevel(zeros(m,n),D,n_level);
dw = U;
bw = U;

%initial guess.
u = nufft_adj(f,st);

dx = zeros(size(u));
dy = zeros(size(u));
bx = zeros(size(u));
by = zeros(size(u));

%movie? sure why not
aviobj = avifile('cs_recon.avi', 'FPS', 10);
aviobj.Quality = 10;
figgn = figure();

%constrained, replace f with fl and update after each unconstrained
fl = f;
errorsr = [0];
errors = [0];
subplot(2,2,1);
imagesc(img);colormap hot;
subplot(2,2,2);

iters = [];

%start the optimization outer loop is constraint enforcement
for ell = 1:100
    %unconstrained
    for k=1:2
        %update u
        rhsD = FraRecMultiLevel(SubFrameletArray(dw,bw),Dt,n_level);
        rhs = mu.*unvec(At(vec(fl)),m,n) + lambda.*Dxt(dx - bx) + lambda.*Dyt(dy - by) + gamma.*rhsD;
        
        %this is where everything sucks.
        [u,flag,reles,iter] = pcg(AtA,vec(rhs),1e-3,2);
        iters = [iters iter];
        
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
        %plot(errors(end - round(length(errors)/2) : end));
        plot(errors);
        subplot(6,2,9);
        plot(errorsr,'r.-');
        %plot(errorsr(end - round(length(errorsr)/2) : end-1),'r.-');

        
        subplot(2,2,4);
        errorfig = abs(u/max(u(:)) - img/max(img(:)));
        imagesc(errorfig);
        
        colormap hot;
        pause(0.01);
    
        fprintf('step ell = %i error (u-img): %f \n', [ell errors(end)]);
        aviobj = addframe(aviobj,figgn);
    end
    
end
close(figgn);
aviobj = close(aviobj);