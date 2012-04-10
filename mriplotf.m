function mriplotf(D)

D1=double(squeeze(D));
DIM = size(D1);
[X,Y,Z]=meshgrid(1:DIM(2),1:DIM(1),1:DIM(3));

slice(X,Y,Z,D1,floor(DIM(1)/2), floor(DIM(2)/2), floor(DIM(3)/2));
colormap(gray);
shading flat;
title('3D Slices')
xlabel('x');ylabel('y');zlabel('z');

end