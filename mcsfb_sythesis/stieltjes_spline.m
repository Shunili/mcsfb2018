function [ xi , pp ] = stieltjes_spline(t,ftt)

n = length(t)-1;
xi    = zeros(n,4);
h = reshape(diff(t)/2, n,1);

pp = pchip(t, ftt);
d = pp.coefs(:,1);
c = pp.coefs(:,2);
e = pp.coefs(:,3);
a = pp.coefs(:,4);
dh3 = d.*h.^3;
ch2 = c.*h.^2;
eh  = e.*h;
xi(:,1) =  5/2*dh3 + 3/2*ch2 + eh + a;
xi(:,2) = 15/4*dh3 +   2*ch2 + eh;
xi(:,3) =  3/2*dh3 + 1/2*ch2;
xi(:,4) =  1/4*dh3;


end

