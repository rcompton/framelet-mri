function d = Deriv(u,dim,t)
% x-direction is derivative along 1st corrd.
% dim specs which direction we diff along
% t is for transpose

shift = zeros(ndims(u),1);

if t
    shift(dim) = -1;
else
    shift(dim) = 1;
end

d = u - circshift(u,shift);

end