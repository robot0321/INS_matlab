%AVP simulator
%attitude, velocity, position

clc;
close all;
clear;

N=1000;
dt = 0.01;
pos = zeros(3,N);
vel = ones(3,N)*1;
acc = zeros(3,N);
gyro = zeros(3,N);
attitude = zeros(3,N);
vel(1,:) = linspace(1,0,N);
vel(3,:) = linspace(0,1,N);

for i=1:N-1
    pos(:,i+1) = pos(:,i) + vel(:,i)*dt;
    acc(:,i) = (vel(:,i+1) - vel(:,i))/dt;
end

gyro(1,:) = ones(1,N)*pi/100;
gyro(3,:) = ones(1,N)*pi/100;
for i=1:N-1
    attitude(:,i+1) = attitude(:,i) + gyro(:,i)*dt;
end


%%

C3 = [cos(attitude(3,i)), sin(attitude(3,i)), 0;
      -sin(attitude(3,i)), cos(attitude(3,i)), 0;
      0, 0, 1];
      
C2 = [cos(attitude(2,i)), 0, -sin(attitude(2,i));
      0, 1, 0;
      sin(attitude(2,i)), cos(attitude(2,i)), 0];
      
C1 = [1, 0, 0;
      0, cos(attitude(1,i)), sin(attitude(1,i));
      0, -sin(attitude(1,i)), cos(attitude(1,i))];
  
  Cb2n = C3*C2*C1;
  Bacc = Cb2n*acc;
  Bgyro = zeros(3,N);
  for i=1:N
     Bgyro(:,i) = [gyro(1,i), 0, 0]' + C3*[0, gyro(2,i), 0]' + C3*C2*[0, 0, gyro(3,i)]'; 
  end
  
  