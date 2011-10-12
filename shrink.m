function [ d ] = shrink( z,tau )

d = sign(z).*max(zeros(size(z)), abs(z) - tau);


end

