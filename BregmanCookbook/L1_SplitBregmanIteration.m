function u=L1_SplitBregmanIteration(f,A,mu,lambda,err)

%======================================
%
% L1 Split Bregman Iteration
% Version:
% - v1.0 - 03/31/2011
%
% u=L1_SplitBregmanIteration(f,A,mu,lambda,err)
%
% This function compute the solution
% u = arg min |u|+0.5*mu||Au-f||_2^2
%
% by using the Split Bregman Iteration
%
% In this version we consider only
% the case where A is a square matrix
%
% f: measured data
% A: some linear operator in its matrix
%    form
% mu: regularization coefficient
% lambda: "spliting" regularization
%         coefficient
%
% Author: Jerome Gilles
% Institution: UCLA - Math Department
% email: jegilles@math.ucla.edu
% 
% Note: typically mu=10, lambda=1,
%       eps=0.01 work well
%======================================
N=size(f,1)

d=zeros(N,1);
b=zeros(N,1);
u=zeros(N,1);

Z=zeros(N,1);
Ft=mu*A'*f;
IV=inv(mu*A'*A+lambda*eye(N));
up=ones(N,1);

while ((u-up)'*(u-up))>eps,
   up=u;
   u=IV*(Ft+lambda*(d-b));
   tmp=u+b;
   d=sign(tmp).*max(Z,abs(tmp)-1/lambda);
   b=tmp-d;
end
