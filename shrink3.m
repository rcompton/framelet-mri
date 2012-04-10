function [xs,ys,zs] = shrink3(x,y,z,lambda)

s = sqrt(x.*conj(x)+y.*conj(y)+z.*conj(z));
ss = s-lambda;
ss = ss.*(ss>0);

s = s+(s<lambda);
ss = ss./s;

xs = ss.*x;
ys = ss.*y;
zs = ss.*z;

end