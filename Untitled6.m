e=0.0818797908425;
Earth_Omega = 7.292115e-5;
Earth_R_short = 6356752.3142;
Earth_R_long = 6378137.0;

z = 0:0.1:Earth_R_short;
x = Earth_R_long*sqrt(1-z.^2/Earth_R_short^2);
zeta = atan(z./x/sqrt(1-e^2));

Lb = atan((z*sqrt(1-e^2) + e^2*Earth_R_long*sin(zeta).^3)/sqrt(1-e^2)./(x-e^2*Earth_R_long*cos(zeta).^3));

plot(Lb)