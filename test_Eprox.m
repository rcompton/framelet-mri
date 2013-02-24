%
% Test matrix decompositions
%
clear all;close all;clc;

%load fmap_test.mat
clear all;close all;clc;

%load fmap_test.mat
fmap = double(flipud(rgb2gray(imread('fmappic.jpg'))));
fmap = fmap/max(fmap(:))*60;

fmap = circshift(fmap,[-10 -4]);
%%
n = 64;

fmap = imresize(fmap,[n n]);

%simulate interpolated field map...
Delta1 = .1;
Arg = @(x) atan(imag(x)./real(x));
ph = phantom('Modified Shepp-Logan',n);

y1 = ph + 1e-6*rand(size(ph));
y2 = ph.*exp(1j*fmap*Delta1) + 1e-6*rand(size(ph));
fmapest = Arg(conj(y1).*y2)/Delta1;

imagesc(fmapest);
colormap bone

%use the estimation
fmap = reshape(fmapest,[n*n 1]);
fmap = reshape(fmap,[n*n 1]);

%for EPI, 2msec/line (ok..)
T = n*(2e-3);
ti = [];
for k=1:n
    ti = [ti linspace(0,T,n)]; %?
end

E = exp(-1j*(fmap*ti));


L = 20;
ntrials = 1;
svderrs = zeros(ntrials,L);
riderrs = zeros(ntrials,L);
for n=1:ntrials
    for l=1:L
    l
    [u,s,v] = pca(E,l,1);
    svderrs(n,l) = norm(E - u*s*v')/norm(E);

    [b,p] = id(E,l);
    riderrs(n,l) = norm(b*p - E)/norm(E);
    end
end

%%
close all;
semilogy(svderrs','b-o','LineWidth',2);
hold on;
semilogy(riderrs','r-x','LineWidth',2);
ylabel('log ||E - E_{approx}|| / ||E||');
xlabel('approximation rank');

legend('SVD','ID');

hAll = findall(gcf);
for idx = 1 : length(hAll)
  try
    set(hAll(idx),'fontsize', 15);
  catch
    % never mind...
  end
end
