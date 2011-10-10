function col=FrameletCoef_to_column(Dec)

L=size(Dec,2); %number of level
P=size(Dec{1},1); %number of coefficient image per level
N=size(Dec{1}{1,1},1); %row number of each image
M=size(Dec{1}{1,1},2); %column number of each image
off=P^2*N*M;

L
P
N
M

L*off

col=zeros(L*off);

for l=1:L
   for px=1:P
      for py=1:P
         for m=1:M
            col([(m-1)*N+1+off:m*N+off])=Dec{l}{px,py}(:,m); 
         end
      end
   end
end