function [B, P] = id(A, k)
%ID Randomized Interpolative Decomposistion
%   Computes a rank-k approximation to A of the form A=B*P where the
%   columns of B are drawn from the columns of A and the entries in P are
%   not too large. A=B*P should be exact on true low-rank matricies.
%   Based on Tygert's PNAS paper.

[m, ~] = size(A);

Y = randn(k,m)*A;

[~,R,PI] = qr(Y,0);
%Q = Q(:,1:k);
R = R(1:k, :);
[~, PIt] = sort(PI);

S = R(1:k,1:k);
T = R(:,k+1:end);
%Z = Q*S;

P = [eye(k) S\T];
P = P(:,PIt);

bestcols = PI(1:k);

B = A(:,bestcols);

end

