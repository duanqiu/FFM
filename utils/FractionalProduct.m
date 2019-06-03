function out = FractionalProduct(V,n)
if (~exist('V','var'))
    V = 0.5;
end
if (~exist('n','var'))
    n = 5;
end

z = zeros((2*n+1),1);

z(n) = V/4 + V^2/8;
z(n+1) = 1- V^2/2 - V^3/8;
z(n+2) = -5*V/4 + 5*V^3/16 + V^4/16;
k = 1;
for i = n+3:2*n-2
    k = k+1;
   z(i) =  (gamma(k-V+1)/gamma(k+2)*(V/4 + V^2/8)+gamma(k-V)/gamma(k+1)*(1- V^2/4)+gamma(k-V-1)/gamma(k)*(-V/4 + V^2/8))/gamma(-V);
end
z(2*n-1) = (gamma(n-V-1)/gamma(n)*(V/4 + V^2/8)+gamma(n-V-2)/gamma(n-1)*(1- V^2/4)+gamma(n-V-3)/gamma(n-2)*(-V/4 + V^2/8))/gamma(-V);
z(2*n) = gamma(n-V-1)/(gamma(n)*gamma(-V))*(1-V^2/4) + gamma(n-V-2)/(gamma(n-1)*gamma(-V))*(-V/4 + V^2/8);
z(2*n+1) = gamma(n-V-1)/(gamma(n)*gamma(-V))*( -V/4 + V^2/8);
out = z;




