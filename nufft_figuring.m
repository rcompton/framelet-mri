%
% nufft example
%
close all;clear all;clc;

addpath(genpath('~./Dropbox/irt'));
addpath('nufft_files/');

%
n1 = 75;
n2 = 100;
npts = 5000
x = imresize(peaks(max(n1,n2)),[n1 n2]);

%the sprial, a 201 point subset of the full freq domain, scale to [0 1] each loc to sprial it
%omega = linspace(0, 20*2*pi, npts)';	% crude spiral:
%omega = pi*[cos(omega) sin(omega)].*omega(:,[1 1])/max(omega);

%all the points...
[k1pts k2pts] = meshgrid(1:n1, 1:n2);
k1pts = reshape(k1pts,numel(k1pts),1)*2*pi/n1;
k2pts = reshape(k2pts,numel(k2pts),1)*2*pi/n2;
omega = [k1pts k2pts];

%random sample of points
omega = omega(sort(randsample(1:max(size(omega)),npts)),:);

%plot the sample locations
plot(omega(:,1), omega(:,2),'b.');

%the number of neighbors to use when interpolating
j1 = 15;
j2 = 15;

%the fft sizes?
k1 = 2*n1;
k2 = 2*n2;

%this sytax
args = {omega, [n1 n2], [j1 j2], [k1 k2]};

%Create fft matrix
%A = Gnufft(args); %fast
%Ad = Gdsft(omega,[n1 n2]); %exact

%old way?
fprintf('nufft init: ');
tic
st = nufft_init(omega, [n1 n2], [j1 j2], [k1 k2]);
toc

%fancy Gnufft way
% tic
% ya = A*x;
% toc

%exact nuft
%tic
%yx = Ad*x;
%toc

%old style?
% tic
% yo = nufft(x,st); %an fft from the whole image into the small, oddly spaced, freq domain
% toc

%scaling? I think nufft is like fftw....
scale_factor = sqrt(npts);
yos = nufft(x,st)/scale_factor;
yosi = nufft_adj(yos,st)/scale_factor;

figure()
subplot(1,2,1)
surf(real(yosi));
subplot(1,2,2)
surf(x);



