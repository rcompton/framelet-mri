function R=AddFrameletArray(A,B)

%=========================================
%
% This function performs the sum of each 
% coefficients between two arrays of
% Framelet coefficients.
% A and B MUST have the same structure
%
% Author: Jerome Gilles
% Institution: UCLA - Math Department
% email: jegilles@math.ucla.edu
%
%=========================================

L=length(A);

for l=1:L
   [NH,NW]=size(A{l});
   for nh=1:NH
      for nw=1:NW
          R{l}{nh,nw}=A{l}{nh,nw}+B{l}{nh,nw};
      end
   end
end