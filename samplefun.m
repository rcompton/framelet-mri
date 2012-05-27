function [ y ] = samplefun(R,B,C,x,t)
%sum Bl*G*Cl x
%R is random sample matrix
%B is n^2 by L
%C is L by n^2
%x is the image, in vector format
%t is for transpose

n = sqrt(size(B,1));
assert(n == round(n));
x = reshape(x,[n n]);
y = zeros(size(x));
scal = sqrt(n*n);

if ~t
    for l=1:size(B,2)
        Bl = reshape(B(:,l),[n n]);
        Cl = reshape(C(l,:),[n n]);
        
        y = y + Bl.*R.*(1/scal).*fft2(Cl.*x);
    end
else
    for l=1:size(B,2)
        Bl = reshape(B(:,l),[n n]);
        Cl = reshape(C(l,:),[n n]);
        
        y = y + conj(Cl).*scal.*ifft2(R.*conj(Bl).*x);
    end
end

y = reshape(y,[numel(y) 1]);
end

