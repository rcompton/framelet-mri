%
% Screw around with giant partial svd
%
%
close all;clear all;

load fmap_test.mat

[m n] = size(fmap);

num_samps = m*n/2.5;

fmap = reshape(fmap,[1 m*n]);
%fmap = 120*fmap./max(fmap(:));%?

T = 1.2;
t = linspace(.01,T,num_samps);


E = zeros(num_samps, m*n);
for i=1:num_samps
    E(i,:) = exp(-fmap*t(i));
end

[B,sigs,C] = pca(E,6);

G = dftmtx(m*n);
G = G(sort(randsample(1:m*n,num_samps)),:);

L = length(diag(sigs));
A = zeros(size(G));
for l=1:L
    A = A + sigs(l,l)*diag(B(:,l))*G*diag(C(:,l));
end


%%
imagesc(imag(G));
caxis([min(imag(G(:))) max(imag(G(:)))]);
colormap jet;
figure()
imagesc(imag(A));
caxis([min(imag(G(:))) max(imag(G(:)))]);
colormap jet;



