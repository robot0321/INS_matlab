clc;
clear;
close all;

M = importdata('data6.txt');
S = zeros(1, size(M,1)-1);

% 1.time // 234 accel // 567 gyro // 8910 mag // 11 tmp
g = 9.8;
ax = M(:,2);
ay = M(:,3);
az = M(:,4);
gx = M(:,5)*pi/180;
gy = M(:,6)*pi/180;
gz = M(:,7)*pi/180;
mx = M(:,8);
my = M(:,9);
mz = M(:,10);
tmp = M(:,11);

abmag = sqrt(ax.^2 + ay.^2 + az.^2);
ax = ax./abmag;
ay = ay./abmag;
az = az./abmag;

% phi = atan(ay./az);
% theta = asin(-ax/1); % g ¥‹¿ß
% ly = mx.*cos(theta)*sin(p%my.*cos(phi) - mz.*sin(phi);
% lz = -mx.*sin(theta) + my.*cos(theta).*sin(phi) + mz.*cos(theta).*cos(phi);
% psi = atan(ly./lz);
phi = atan(ay./sqrt(ax.^2+az.^2));
theta = atan(ax./sqrt(ay.^2+az.^2));
psi = atan(az./sqrt(ax.^2+az.^2));


for i=2:size(phi)
   if(abs(phi(i)+pi - phi(i-1)) < abs(phi(i)-phi(i-1)))
      phi(i) = phi(i) + pi;
   end
   if(abs(phi(i)-pi - phi(i-1)) < abs(phi(i)-phi(i-1)))
       phi(i) = phi(i) - pi;
   end
end

for i=2:size(psi)
   if(abs(psi(i)+pi - psi(i-1)) < abs(psi(i)-psi(i-1)))
      psi(i) = psi(i) + pi;
   end
   if(abs(psi(i)-pi - psi(i-1)) < abs(psi(i)-psi(i-1)))
       psi(i) = psi(i) - pi;
   end
end


figure(1);
subplot(3,1,1);
plot(phi)
title('phi');
subplot(3,1,2);
plot(theta);
title('theta');
subplot(3,1,3);
plot(psi);
title('psi');


dt = 0.01;
grpy = zeros(3, size(M,1));

DCM = zeros(3,3,size(M,1));
DCM(:,:,1) = eye(3);
EE = zeros(3,size(M,1));
for i=1:size(M,1)-1
    %omega = [1, -gz(i)*dt, gy(i)*dt; gz(i)*dt, 1, -gx(i)*dt; -gy(i)*dt, gx(i)*dt, 1];
    %DCM(:,:,i+1) = DCM(:,:,i)*omega;
    %rot = [1, sin(grpy(1,i))*tan(grpy(2,i)), cos(grpy(1,i))*tan(grpy(2,i)); 0, cos(grpy(1,i)), -sin(grpy(1,i)); 0, sin(grpy(1,i))/cos(grpy(2,i)), cos(grpy(1,i))/cos(grpy(2,i))];
    grpy(:,i+1) = grpy(:,i) + [(gy(i)*sin(grpy(1,i))+gz(i)*cos(grpy(1,i)))*tan(grpy(2,i))+gx(i);gy(i)*cos(grpy(1,i))-gz(i)*sin(grpy(1,i));(gy(i)*sin(grpy(1,i))+gz(i)*cos(grpy(1,i)))*sec(grpy(2,i))]*dt;%+ rot*[gx(i);gy(i);gz(i)];
    %[yy,pp,rr] = dcm2angle(DCM(:,:,i+1));
    %EE(:,i+1) = [DCM(2,3,i+1); DCM(3,1,i+1); DCM(1,2,i+1)];%[yy,pp,-rr];
%     aa = arccos((DCM(1,1,i+1)+DCM(2,2,i+1)+DCM(3,3,i+1)-1)/2);
%     roll = (DCM(3,2,i+1) - DCM(2,3,i+1))/2/sin(aa);
%     pitch = (DCM(1,3,i+1) - DCM(3,1,i+1))/2/sin(aa);
%     yaw = (DCM(2,1,i+1) - DCM(1,2,i+1))/2/sin(aa);
end

% gphi = zeros(1,size(M,1))';
% gpsi = zeros(1,size(M,1))';
% gtheta = zeros(1,size(M,1))';
% % drpy = [1, sin(phi).*tan(theta), cos(phi).*tan(theta); 0, cos(phi), -sin(phi); 0, sin(phi)./cos(theta), cos(phi)./cos(theta)]*[M(:,5)'; M(:,6)'; M(:,7)'];
% % grpy = grpy + drpy*dt;
% 
% dphi = M(:,5) + sin(gphi).*tan(gtheta).*M(:,6) + cos(gphi).*tan(gtheta).*M(:,7);
% dtheta = M(:,6).*cos(gphi) - M(:,7).*sin(gphi);
% dpsi = M(:,6).*sin(gphi)./cos(gtheta) + M(:,7).*cos(gphi)./cos(gtheta);
% 
% for i=1:size(M,1)-1
%     dphi = M(i,5) + sin(gphi(i)).*tan(gtheta(i)).*M(i,6) + cos(gphi(i)).*tan(gtheta(i)).*M(i,7);
%     dtheta = M(i,6).*cos(gphi(i)) - M(i,7).*sin(gphi(i));
%     dpsi = M(i,6).*sin(gphi(i))./cos(gtheta(i)) + M(i,7).*cos(gphi(i))./cos(gtheta(i));
%     gphi(i+1) = gphi(i) + dphi*dt;
%     gpsi(i+1) = gpsi(i) + dpsi*dt;
%     gtheta(i+1) = gtheta(i) + dtheta*dt;
% end



figure(2);
subplot(3,1,1);
plot(grpy(1,:));title('phi');
subplot(3,1,2);
plot(grpy(2,:));title('theta');
subplot(3,1,3);
plot(grpy(3,:));title('psi');

%%
figure(3);
a = [zeros(size(M,1),1), theta, phi];
%yaw_c = 0;
yaw2_c = psi(1);
pitch_c = theta(1);
roll_c = phi(1);
yaw_g = grpy(3,1);%gpsi(1);
pitch_g = grpy(2,1);%gtheta(1);
roll_g = grpy(1,1);%gphi(1);

for frame=1:1:size(M,1)
    %Read Data;
    %Convert Quaternion to Roll Pitch Yaw
    yaw = 0;
    yaw2 = psi(frame);
    pitch = theta(frame);
    roll = phi(frame);
    yaw_gc = grpy(3,frame);%gpsi(frame);
    pitch_gc = grpy(2,frame);%gtheta(frame);
    roll_gc = grpy(1,frame);%gphi(frame); 
    
    %Visualize Data On Sphere Or any Other Objects
    [x,y,z] = sphere;
 
    subplot(1,3,1);
    h = surf(x,y,z);axis('square'); 
    title('only acc')
    xlabel('x'); ylabel('y'); zlabel('z');
    %Rotate Object
    rotate(h,[1,0,0],(roll-roll_c)*180/pi)
    rotate(h,[0,1,0],(pitch-pitch_c)*180/pi)
    %rotate(h,[0,0,1],(yaw-yaw_c)*180/pi+90)
    rotate(h,[0,0,1],(0)*180/pi)
    
    subplot(1,3,2);
    h = surf(x,y,z);axis('square'); 
    title('acc mag')
    xlabel('x'); ylabel('y'); zlabel('z');
    %Rotate Object
    rotate(h,[1,0,0],(roll-roll_c)*180/pi)
    rotate(h,[0,1,0],(pitch-pitch_c)*180/pi)
    %rotate(h,[0,0,1],(yaw-yaw_c)*180/pi)
    rotate(h,[0,0,1],(yaw2-yaw2_c)*180/pi+90)
    %view(0,0);
    
    subplot(1,3,3);
    h = surf(x,y,z);axis('square'); 
    title('gyro')
    xlabel('x'); ylabel('y'); zlabel('z');
    %Rotate Object
    rotate(h,[0,0,1],(yaw_gc - yaw_g)*180/pi)
    
    rotate(h,[0,1,0],(pitch_gc - pitch_g)*180/pi)
    rotate(h,[1,0,0],(roll_gc - roll_g)*180/pi)
    %rotate(h,[0,0,1],(yaw-yaw_c)*180/pi+90)
    
    %view(0,0);
    
    drawnow
end