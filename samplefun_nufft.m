function [ y ] = samplefun_nufft(st,B,C,x,m,n,t)
%sum Bl*G*Cl x
%x is a column vector (you need this for pcg etc.)
%output is a column vector
%t is for transpose
%For the forward nufft, the input is a 2D image which outputs to a 1D
%array with num_samples elements. For the adjoint nufft the input is 1D and
%the output is 2D
%
%C is a cell array of L matrices that you use for the inhomogenity
%correction. B is the same deal
%Everything column vectors! Like in the papers.

[mn,~] = size(C{1});
num_samples = length(B{1});


if ~t
    assert(mn ==  length(x));
else
    assert(num_samples == length(x));
end
assert(n == mn/m);

L = length(C);
assert(length(B) == L);%gawd ryan who cares?

scale = sqrt(m*n);

if t
    y = zeros(m,n);
    for l=1:L
        Bx = B{l}.*x;
        y = y + reshape(C{l},[m n]).*nufft_adj(Bx, st)/scale;
    end
    y = reshape(y, [numel(y) 1]);
    
else
    y = zeros(num_samples,1);

    for l=1:L
        Cx = C{l}.*x;
        y = y + B{l}.*nufft(reshape(Cx,[m n]), st)/scale;
    end
end


end
