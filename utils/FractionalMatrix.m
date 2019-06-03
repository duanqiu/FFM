function [MY,MX] = FractionalMatrix(N,v1,r)


if (~exist('r','var'))
    r = 10;
end
if (~exist('N','var'))
    N = 110;
end
if (~exist('v1','var'))
    v1 = 0.5;
end


VV1 = FractionalProduct(v1);

ope_size = length(VV1);
% ope_size = length(11);
B = VV1;
rr = floor(ope_size/2);


MY = spdiags(zeros(N,1),0,N,N);
for i=1:ope_size
	MY = MY + spdiags(B(i).*ones(N,1),i-rr-1,N,N);
end

MX = spdiags(zeros(N,1),0,N,N);
for i=1:ope_size
	MX = MX + spdiags(B(i).*ones(N,1),(i-rr-1)*r,N,N);
end