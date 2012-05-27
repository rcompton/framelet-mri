%
% Directly make a 64^2x64^2 E matrix
%

clear all;close all;clc;

load fmap_test.mat

n = 128;

fmap = imresize(fmap,[n n]);
fmap = rot90(fmap,3);

%simulate interpolated field map...
Delta1 = .1;
Arg = @(x) atan(imag(x)./real(x));
ph = phantom('Modified Shepp-Logan',n);

y1 = ph + 1e-1*rand(size(ph));
y2 = ph.*exp(1j*fmap*Delta1) + 1e-1*rand(size(ph));
fmapest = Arg(conj(y1).*y2)/Delta1;

imagesc(fmapest);
colormap bone

%use the estimation
fmap = reshape(fmapest,[n*n 1]);
fmap = reshape(fmap,[n*n 1]);

%for EPI, 2msec/line (ok..)
T = n*(10e-2);
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
save onetwentyeightE;

%%
test_fmap_blur