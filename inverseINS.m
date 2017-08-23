%inverse INS (local)
clc;
clear;
close all;

dt = 0.01; %100Hz
Earth_Omega = 7.292115e-5;
Earth_R_short = 6356752.3142;
Earth_R_long = 6378137.0;
flattening_f = (Earth_R_long - Earth_R_short)/Earth_R_long;
eccentricity_e = sqrt(flattening_f*(2-flattening_f));
GM = 398600.4418*1000^3;
[ex, ey, ez] = ellipsoid(0,0,0,Earth_R_long,Earth_R_long,Earth_R_short,20);

% %직진
% 위치: 위도(L), 경도(l), 중심부터거리(r)
% init_P = [36.0*pi/180, 127.0*pi/180, Earth_R_long];
% final_P = [37.0*pi/180, 127.0*pi/180, Earth_R_long];
% init_C_n_b = [1, 0, 0; 0, 1, 0; 0, 0, 1];
% P1 = init_P(1):0.0001:final_P(1);
% P1(2,:) = ones(1,size(P1,2))*init_P(2);
% P1(3,:) = ones(1,size(P1,2))*init_P(3);
% C_n_b = zeros(3,3,size(P1,2));
% for i=1:size(P1,2)
%     C_n_b(:,:,i) = [1, 0, 0; 0, 1, 0; 0, 0, 1];
% end
% r_mag = R_surface(Earth_R_long, Earth_R_short, P1(1,:));

% %회전 
% wr = 0.01; % 회전 반경 (rad)
% init_P = [36.0*pi/180 + wr*cos(0), 127.0*pi/180 + wr*sin(0), Earth_R_long];
% init_C_n_b = [1, 0, 0; 0, 1, 0; 0, 0, 1];
% Ns = 10000;
% radi = linspace(0, 2*pi, Ns);
% P1 = 36.0*pi/180 + wr*cos(radi);
% P1(2,:) = 127.0*pi/180 + wr*sin(radi);
% P1(3,:) = init_P(3) + 500*sin(radi);%ones(1,Ns) * init_P(3);
% C_n_b = zeros(3,3,size(P1,2));
% for i=1:size(P1,2)
%     C_n_b(:,:,i) = [1, 0, 0; 0, 1, 0; 0, 0, 1];
% end
% r_mag = R_surface(Earth_R_long, Earth_R_short, P1(1,:));%*(1+tan(wr)^2);


%정지 (동체 회전)
%위치: 위도(L), 경도(l), 중심부터거리(r)
init_P = [36.0*pi/180, 127.0*pi/180, Earth_R_long];
init_C_n_b = [1, 0, 0; 0, 1, 0; 0, 0, 1];
Ns = 100000;
P1 = ones(1,Ns)*init_P(1);
P1(2,:) = ones(1,Ns)*init_P(2);
P1(3,:) = ones(1,Ns)*init_P(3);
C_n_b(:,:,1) = init_C_n_b;
wb = [0.3, 0.3, 0.3]';
omega = [0, -wb(3), wb(2); wb(3), 0, -wb(1); -wb(2), wb(1), 0];
for i=1:size(P1,2)-1
    C_n_b(:,:,i+1) = C_n_b(:,:,i)*(eye(3) + omega*dt);
end



% %직진 (동체 회전)
% %위치: 위도(L), 경도(l), 중심부터거리(r)
% init_P = [36.0*pi/180, 127.0*pi/180, Earth_R_long];
% final_P = [37.0*pi/180, 127.0*pi/180, Earth_R_long];
% init_C_n_b = [1, 0, 0; 0, 1, 0; 0, 0, 1];
% P1 = init_P(1):0.0001:final_P(1);
% P1(2,:) = ones(1,size(P1,2))*init_P(2);
% P1(3,:) = ones(1,size(P1,2))*init_P(3);
% C_n_b(:,:,1) = init_C_n_b;
% wb = [0.3, 0.3, 0.3]';
% omega = [0, -wb(3), wb(2); wb(3), 0, -wb(1); -wb(2), wb(1), 0];
% for i=1:size(P1,2)-1
%     C_n_b(:,:,i+1) = C_n_b(:,:,i)*(eye(3) + omega*dt);
% end
% r_mag = R_surface(Earth_R_long, Earth_R_short, P1(1,:));


% %회전 (동체 회전)
% wr = 0.01; % 회전 반경 (rad)
% init_P = [36.0*pi/180 + wr*cos(0), 127.0*pi/180 + wr*sin(0), Earth_R_long];
% init_C_n_b = [1, 0, 0; 0, 1, 0; 0, 0, 1];
% Ns = 10000;
% radi = linspace(0, 2*pi, Ns);
% P1 = 36.0*pi/180 + wr*cos(radi);
% P1(2,:) = 127.0*pi/180 + wr*sin(radi);
% P1(3,:) = init_P(3) + 5000*sin(radi);%ones(1,Ns) * init_P(3);
% C_n_b(:,:,1) = init_C_n_b;
% wb = [1, 1, 1]'*1;
% omega = [0, -wb(3), wb(2); wb(3), 0, -wb(1); -wb(2), wb(1), 0];
% for i=1:size(P1,2)-1
%     C_n_b(:,:,i+1) = C_n_b(:,:,i)*(eye(3) + omega*dt);
%     C_n_b(:,:,i+1) = DCM_ortho_normal_compensation(C_n_b(:,:,i+1));%DCM orthogonalization, normalization
% end
% r_mag = R_surface(Earth_R_long, Earth_R_short, P1(1,:));%*(1+tan(wr)^2);


%% 역 연산
%R_surface(Earth_R_long, Earth_R_short, P1(1,:));%sqrt((1+tan(P1(1,:)).^2)./(Earth_R_short^2 + (Earth_R_long^2)*tan(P1(1,:)).^2))*Earth_R_short*Earth_R_long + P1(3,i);
%R_vector = r_mag.*[cos(P1(1,:)).*cos(P1(2,:)); cos(P1(1,:)).*sin(P1(2,:)); sin(P1(1,:))]; %P1(3,:).*
R_vector = [P1(3,1:end).*cos(P1(1,1:end)).*cos(P1(2,1:end)); P1(3,1:end).*cos(P1(1,1:end)).*sin(P1(2,1:end)); P1(3,1:end).*sin(P1(1,1:end))];
%g = GM*1./r_mag.^2;
g0 = 9.780318 * (1 + 5.3024e-3 .* sin(P1(1,:)).*sin(P1(1,:)) - 5.9e-6 .* sin(2*P1(1,:)) .* sin(2*P1(1,:)));
g = g0./(1+(P1(3,:)-Earth_R_long)/Earth_R_long);
w_n_ie = [Earth_Omega*cos(P1(1,:)); zeros(size(P1(1,:))); -Earth_Omega*sin(P1(1,:))];
g_n_l = [zeros(size(g)); zeros(size(g)); g] - cross(w_n_ie, cross(w_n_ie, R_vector)); %애초에 중력을 일케 쓰면 무조건 로컬

delta_L=zeros(1,size(P1,2));delta_l=zeros(1,size(P1,2));delta_h=zeros(1,size(P1,2));
for i=1:size(P1,2)-1
    delta_L(i) = (P1(1,i+1) - P1(1,i))/dt; 
    delta_l(i) = (P1(2,i+1) - P1(2,i))/dt; 
    delta_h(i) = (P1(3,i+1) - P1(3,i))/dt; 
end


w_n_en = [delta_l.*cos(P1(1,:)); -delta_L; -delta_l.*sin(P1(1,:))];
w_n_in = w_n_ie + w_n_en;
w_b_in = zeros(3, size(P1,2));
for i=1:size(P1,2)
    w_b_in(:,i) = C_n_b(:,:,i)' * w_n_in(:,i);
end

delta_C_n_b = zeros(3,3,size(P1,2));
Omega_b_nb = zeros(3,3,size(P1,2));
w_b_nb = zeros(3, size(P1,2));
for i=1:size(P1,2)-1
    delta_C_n_b(:,:,i) = (C_n_b(:,:,i+1) - C_n_b(:,:,i))/dt;
end
for i=1:size(P1,2)
    Omega_b_nb(:,:,i) = C_n_b(:,:,i)'*delta_C_n_b(:,:,i);
    w_b_nb(:,i) = wb;%[Omega_b_nb(3,2,i); Omega_b_nb(1,3,i); Omega_b_nb(2,1,i)];   
end
w_b_ib = w_b_in + w_b_nb; %자이로 센서 데이터

V_N = delta_L.*(P1(3,:));
V_E = delta_l.*(P1(3,:)).*cos(P1(1,:));
V_D = -delta_h;
V_n = [V_N; V_E; V_D];
init_V_n = V_n(:,1);

delta_V_n = zeros(3, size(P1,2));
for i=1:size(P1,2)-1
   delta_V_n(:,i) = (V_n(:,i+1)-V_n(:,i))/dt;
end

f_n = delta_V_n + cross(2*w_n_ie + w_n_en, V_n) - g_n_l;
f_b = zeros(3, size(P1,2));
for i=1:size(P1,2)
   f_b(:,i) = C_n_b(:,:,i)'*f_n(:,i);
end

true_V = V_n;
true_P = P1;
true_C_n_b = C_n_b;

save body_data.mat f_b w_b_ib init_P init_C_n_b init_V_n true_P true_V true_C_n_b;

figure(1);
subplot(3,1,1);
plot(P1(1,1:end-1)*180/pi,'LineWidth', 2)
subplot(3,1,2);
plot(P1(2,1:end-1)*180/pi,'LineWidth', 2)
subplot(3,1,3);
plot(P1(3,:),'LineWidth', 2)
figure(2);
% surf(ex, ey, ez);
% axis equal
% hold on;
plot3(P1(3,1:end-1).*cos(P1(1,1:end-1)).*cos(P1(2,1:end-1)), P1(3,1:end-1).*cos(P1(1,1:end-1)).*sin(P1(2,1:end-1)), P1(3,1:end-1).*sin(P1(1,1:end-1)), 'LineWidth', 1, 'Marker', '+')
figure(3);
plot(P1(1,:)*180/pi, P1(2,:)*180/pi,'LineWidth', 2);


% M = C_n_b;
% for frame=1:1:size(M,3)
%     %Read Data;
%     %Convert Quaternion to Roll Pitch Yaw
%     yaw2 = psi(frame);
%     pitch = theta(frame);
%     roll = phi(frame);
%     
%     %Visualize Data On Sphere Or any Other Objects
%     [x,y,z] = sphere;
%  
%     subplot(1,3,1);
%     h = surf(x,y,z);axis('square'); 
%     title('only acc')
%     xlabel('x'); ylabel('y'); zlabel('z');
%     %Rotate Object
%     rotate(h,[1,0,0],(roll-roll_c)*180/pi)
%     rotate(h,[0,1,0],(pitch-pitch_c)*180/pi)
%     %rotate(h,[0,0,1],(yaw-yaw_c)*180/pi+90)
%     rotate(h,[0,0,1],(0)*180/pi)
%     
%     subplot(1,3,2);
%     h = surf(x,y,z);axis('square'); 
%     title('acc mag')
%     xlabel('x'); ylabel('y'); zlabel('z');
%     %Rotate Object
%     rotate(h,[1,0,0],(roll-roll_c)*180/pi)
%     rotate(h,[0,1,0],(pitch-pitch_c)*180/pi)
%     %rotate(h,[0,0,1],(yaw-yaw_c)*180/pi)
%     rotate(h,[0,0,1],(yaw2-yaw2_c)*180/pi+90)
%     %view(0,0);
%     
%     subplot(1,3,3);
%     h = surf(x,y,z);axis('square'); 
%     title('gyro')
%     xlabel('x'); ylabel('y'); zlabel('z');
%     %Rotate Object
%     rotate(h,[0,0,1],(yaw_gc - yaw_g)*180/pi)
%     
%     rotate(h,[0,1,0],(pitch_gc - pitch_g)*180/pi)
%     rotate(h,[1,0,0],(roll_gc - roll_g)*180/pi)
%     %rotate(h,[0,0,1],(yaw-yaw_c)*180/pi+90)
%     
%     %view(0,0);
%     
%     drawnow
% end

