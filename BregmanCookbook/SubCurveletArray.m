function R=SubCurveletArray(A,B)

%=========================================
%
% This function performs the difference of
% each Curvelet coefficients between two
% arrays of Curvelet coefficients.
% A and B MUST have the same structure
%
% Author: Jerome Gilles
% Institution: UCLA - Math Department
% email: jegilles@math.ucla.edu
%
%=========================================

L=length(A);

for l=1:L
   N=length(A{l});
   for nh=1:N
          R{l}{nh}=A{l}{nh}-B{l}{nh};
   end
end