clc;
clear;
close all;
figure();
axis([-3.16e6, -3.06e6, 4e6, 4.15e6, 3.74e6, 3.84e6]);

dt = 0.01; %100Hz
Earth_Omega = 7.292115e-5;
Earth_R_short = 6356752.3142;
Earth_R_long = 6378137.0;
flattening_f = (Earth_R_long - Earth_R_short)/Earth_R_long;
eccentricity_e = sqrt(flattening_f*(2-flattening_f));
[ex, ey, ez] = ellipsoid(0,0,0,Earth_R_long,Earth_R_long,Earth_R_short,20);

R = 7000000; %(m)
Pe = [36.0, 127.0;
       37.0, 127.0;
       36.5, 128.0;
       36.0, 127.0];
h = R_surface(Earth_R_long, Earth_R_short, Pe(1,1));

for i=1:size(Pe,1)-1
   Lati = linspace(Pe(i,1), Pe(i+1,1), 100)';
   Long = linspace(Pe(i,2), Pe(i+1,2), 100)';
   for j=1:100
       heig(j,1) = R_surface(Earth_R_long, Earth_R_short, Lati(j)*pi/180);
   end
   Pc(100*i-99:100*i, 1:3) = [heig.*cos(Lati*pi/180).*cos(Long*pi/180), heig.*cos(Lati*pi/180).*sin(Long.*pi/180), heig.*sin(Lati*pi/180)];
end

%surf(ex,ey,ez);
hold on;
plot3(Pc(1:end,1),Pc(1:end,2),Pc(1:end,3));
pause(1);
hold off;

figure();
for i=1:300
plot3(Pc(1:i,1),Pc(1:i,2),Pc(1:i,3),'Color','red','LineWidth', 1);
axis([-3.16e6, -3.06e6, 4e6, 4.15e6, 3.74e6, 3.84e6]);
pause(0.001);
end