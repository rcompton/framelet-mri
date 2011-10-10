% Let img be the image. First use
%
% [D,R]=GenerateFrameletFilter(1)
%
% to generate the tight frame filters. 
% 
% Then use
%
% coeff=FraDecMultiLevel(img,D,Level)
%
% to decompose.
%
% To thresh the coefficient, you should first use
%
% muLevel=getwThresh(mu,-1,Level,D);
%
% to generate thresholdings for different levels; then use
%
% tcoeff=CoeffOper('s',coeff,muLevel);
%
% You can also use CoeffOper to perform operatations to the coefficients such as + - 
%
% Finally, the reconstruction can be done by
%
% rimg=FraRecMultiLevel(tcoeff,R,Level)