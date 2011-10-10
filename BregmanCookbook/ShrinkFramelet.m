function R=ShrinkFramelet(A,tau)

%===================================================
%
%  Execute the shrinkage of Framelet coefficients
%  stored in array A with threshhold tau. The low
%  "frequencies" coefficients are keep without any
%  changes.
%  
%  Parameters:
%       A: input array of Framelet coeeficients
%       tau: threshhold
%       R: output array of threshholded coefficients
%
%   Author: Jerome Gilles
%   Institution: UCLA - Math Department
%   email: jegilles@math.ucla.edu
%====================================================

L=length(A);

for l=1:L
   [NH,NW]=size(A{l});
   for nh=1:NH
      for nw=1:NW
          if ((nw == 1) && (nh == 1))
             R{l}{nh,nw}=A{l}{nh,nw};
          else
             R{l}{nh,nw}=sign(A{l}{nh,nw}).*max(zeros(size(A{l}{nh,nw})),abs(A{l}{nh,nw})-tau);
          end
      end
   end
end