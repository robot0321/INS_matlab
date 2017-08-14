%data
close all;
clear all;

dt=0.1;
t=0:dt:75-dt;
T = 5;
N=length(t);
cir = 30*sin(2*pi/T*t);

roll = zeros(1,N*T);
pitch = zeros(1,N*T);
yaw = zeros(1,N*T);
roll(N+1:2*N) = cir;
roll(3*N+1:4*N) = cir;
pitch(2*N+1:3*N) = cir;

noise = 0.05 + 0.01*randn(3,N*T);

roll = roll + noise(1,:);
pitch = pitch + noise(2,:);
yaw = yaw + noise(3,:);

figure();
tt=0:dt:T*75-dt;
plot(tt, roll,'LineWidth',3.0);
hold on;
plot(tt, pitch,'LineWidth',2.0);
plot(tt, yaw);
hold off;

%%

esave = zeros(N*T,3);
for i=1:N*T
    esave(i,:) = EulerGyro(roll(i)*pi/180, pitch(i)*pi/180, yaw(i)*pi/180, dt)*180/pi;
end

figure();
hold on;
subplot(3,1,1);
plot(tt, esave(:,1),'LineWidth',3.0);
subplot(3,1,2);
plot(tt, esave(:,2),'LineWidth',2.0);
subplot(3,1,3);
plot(tt, esave(:,3),'LineWidth',1.0);
hold off;