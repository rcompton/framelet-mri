%
% Test out what the field map does to my image
%
clear all; close all; clc;

load mri;
load onetwentyeightE.mat

img = double(D(:,:,10));
img = imresize(img,[n n]);

img = phantom('Modified Shepp-Logan',n);

%EF = (B*C).*kron(dftmtx(n),dftmtx(n));
%imgblr = ifft2(reshape( EF*reshape(img,[n*n 1]), [n n]) );

%EFimg = EF*reshape(img,[n*n 1]);

%F = kron(dftmtx(n),dftmtx(n)); %2D dft matrix for stacked arrays...

BCimg = fftshift(img);
for l=1:size(B,2)
    Bl = reshape(B(:,l),[n n]);
    Cl = reshape(C(l,:),[n n]);
    
    BCimg = BCimg + Bl.*fft2(Cl.*img);
    
end

imgfastblr = ifft2(BCimg);

figure
imagesc(img);
colormap bone
%figure
%imagesc(abs(imgblr));
figure
imagesc(abs(imgfastblr));
colormap bone