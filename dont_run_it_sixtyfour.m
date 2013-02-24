%
% Directly make a 64^2x64^2 E matrix
%

clear all;close all;clc;

%load fmap_test.mat
fmap = double(flipud(rgb2gray(imread('fmappic.jpg'))));
fmap = fmap/max(fmap(:))*60;

fmap = circshift(fmap,[-10 -4]);

n = 128;

fmap = imresize(fmap,[n n]);

%simulate interpolated field map...
Delta1 = .1;
Arg = @(x) atan(imag(x)./real(x));
ph = phantom('Modified Shepp-Logan',n);


load mri;
ph = double(D(:,:,19));
ph = ph/norm(ph);

%load raw' data.mat';
%ph = flipud(fftshift(abs(fft2(brain))));
%ph = imresize(ph,[128 128]);
%ph = ph/norm(ph);
%ph = ph.*(ph>.003);

y1 = ph + 0e-3*randn(size(ph));
y2 = ph.*exp(1j*fmap*Delta1) + 0e-3*randn(size(ph));
fmapest = Arg(conj(y1).*y2)/Delta1;


addpath('./cm_and_cb_utilities/');
hh = figure();
imagesc(fmapest);
colormap bone
axis off
%colorbar
%cblabel('Hz');

hAll = findall(gcf);
for idx = 1 : length(hAll)
  try
    set(hAll(idx),'fontsize', 15);
  catch
    % never mind...
  end
end

%% 3D fmap
close all;

load twentytwomri128.mat
ph3 = deenoise;

fmap_3d = zeros(n,n,n);
for k=1:n
    fmap_3d(:,:,k) = fmapest*(1 - k/n)^2*(k/n)^2;
end

Delta13 = .3;
y1 = ph3 + 4e-3*randn(size(ph3));
y2 = ph3.*exp(1j*fmap_3d*Delta13) + 4e-3*randn(size(ph3));
fmapest_3d = Arg(conj(y1).*y2)/Delta13;


whitebg('k');
h = vol3d('cdata',fmapest_3d.*ph3,'texture','3D');

view(3);
axis tight;
%daspect([1 1 .4])
alphamap('rampup');
colormap bone;
%%
alphamap(.5 .* alphamap);

%%
%use the estimation
fmap = reshape(fmapest,[n*n 1]);
fmap = reshape(fmap,[n*n 1]);

%for EPI, 2msec/line (ok..)
T = n*(2e-4);
ti = [];
for k=1:n
    ti = [ti linspace(0,T,n)]; %?
end


%%
%E_{ij} = exp(-i(w_j t_i))
E = exp(-1j*(fmap*ti));


[B, C] = id(E, 10);

normest(E-B*C)/normest(E)

clear E;
save('onetwentyeightE','B','C');

%%
%test_fmap_blur