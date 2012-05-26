%
% Test matrix decompositions
%
clear all;close all;clc;

load fmap_test.mat

fmap = imresize(fmap,[64 64]);

fmap = single(fmap);
fmap = fmap + 1j*fmap.*rand(size(fmap));


ts = linspace(0,1,numel(fmap));
fmap = reshape(fmap,numel(fmap),1);

E = exp(-fmap*ts);

L = 20;
ntrials = 1;
svderrs = zeros(ntrials,L);
riderrs = zeros(ntrials,L);
for n=1:ntrials
    for l=1:L

    [u,s,v] = pca(E,l,1);
    svderrs(n,l) = norm(E - u*s*v')/norm(E);

    [b,p] = id(E,l);
    riderrs(n,l) = norm(b*p - E)/norm(E);
    end
end


semilogy(svderrs','b-o');
hold on;
semilogy(riderrs','r-x');