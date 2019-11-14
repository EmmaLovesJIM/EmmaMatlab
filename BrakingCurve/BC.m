s = linspace(1,1000);
v0 = 44;
a = -1*v0^2./max(2*s);

v = sqrt(2*a*s+v0^2);

plot(s,v)