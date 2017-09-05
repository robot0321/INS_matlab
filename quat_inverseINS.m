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
e = 0.0818191908425;

%정지 (동체 회전)
init_P = [36.0*pi/180, 127.0*pi/180, 0];
init_C_n_b = [1, 0, 0; 0, 1, 0; 0, 0, 1];
init_q = [1;0;0;0];
Ns = 20000;
P1 = init_P(1)*ones(1,Ns);
P1(2,:) = init_P(2)*ones(1,Ns);
P1(3,:) = init_P(3)*ones(1,Ns);
C_n_b(:,:,1) = double(init_C_n_b);
qC(:,:,1) = double(init_C_n_b);
q(:,1) = double(init_q);
w_b_nb = [0, 0, 10]'*ones(1,size(P1,2));
%omega = [0, -w_nb(3), w_nb(2); w_nb(3), 0, -w_nb(1); -w_nb(2), w_nb(1), 0];
for i=1:size(P1,2)-1
    %C_n_b(:,:,i+1) = C_n_b(:,:,i)*(eye(3) + omega*dt);
    %C_n_b(:,:,i+1)= DCM_ortho_normal_compensation(C_n_b(:,:,i+1));%DCM orthogonalization, normalization
     q(:,i+1) = quat_update(q(:,i),w_b_nb(:,i),dt);
     q(:,i+1) = q(:,i+1)/norm(q(:,i+1)); %normalization
     qC(:,:,i+1) = quat2DCM(q(:,i+1)); %DCM
    %qC(:,:,i+1) = init_C_n_b;
end

r_mag = R_surface(Earth_R_long, Earth_R_short, P1(1,:));%*(1+tan(wr)^2);


% %직진 (동체 회전)
% %위치: 위도(L), 경도(l), 중심부터거리(r)
% init_P = [36.0*pi/180, 127.0*pi/180, Earth_R_long];
% final_P = [37.0*pi/180, 127.0*pi/180, Earth_R_long];
% init_C_n_b = [1, 0, 0; 0, 1, 0; 0, 0, 1];
% init_q = [1;0;0;0];
% P1 = init_P(1):0.0001:final_P(1);
% P1(2,:) = ones(1,size(P1,2))*init_P(2);
% P1(3,:) = ones(1,size(P1,2))*init_P(3);
% C_n_b(:,:,1) = init_C_n_b;
% qC(:,:,1) = init_C_n_b;
% q(:,1)=init_q;
% w_b_nb = [0.3, 0.3, 0.3]'*ones(1,size(P1,2));
% for i=1:size(P1,2)-1
%      q(:,i+1) = quat_update(q(:,i),w_b_nb(:,i),dt);
%      q(:,i+1) = q(:,i+1)/norm(q(:,i+1)); %normalization
%      qC(:,:,i+1) = quat2DCM(q(:,i+1)); %DCM
% end
% r_mag = R_surface(Earth_R_long, Earth_R_short, P1(1,:));


% %회전 (동체 회전)
% wr = 0.01; % 회전 반경 (rad)
% init_P = [36.0*pi/180 + wr*cos(0), 127.0*pi/180 + wr*sin(0), Earth_R_long];
% init_C_n_b = [1, 0, 0; 0, 1, 0; 0, 0, 1];
% init_q = [1;0;0;0];
% Ns = 10000;
% radi = linspace(0, 8*pi, Ns);
% P1 = 36.0*pi/180 + wr*cos(radi);
% P1(2,:) = 127.0*pi/180 + wr*sin(radi);
% P1(3,:) = init_P(3) + 5000*sin(radi);%ones(1,Ns) * init_P(3);
% C_n_b(:,:,1) = init_C_n_b;
% qC(:,:,1) = init_C_n_b;
% q(:,1)=init_q;
% w_b_nb = [0, 0, 10]'*ones(1,size(P1,2));
% %omega = [0, -w_nb(3), w_nb(2); w_nb(3), 0, -w_nb(1); -w_nb(2), w_nb(1), 0];
% for i=1:size(P1,2)-1
%     %C_n_b(:,:,i+1) = C_n_b(:,:,i)*(eye(3) + omega*dt);
%     %C_n_b(:,:,i+1)= DCM_ortho_normal_compensation(C_n_b(:,:,i+1));%DCM orthogonalization, normalization
%      q(:,i+1) = quat_update(q(:,i),w_b_nb(:,i),dt);
%      q(:,i+1) = q(:,i+1)/norm(q(:,i+1)); %normalization
%      qC(:,:,i+1) = quat2DCM(q(:,i+1)); %DCM
%     %qC(:,:,i+1) = init_C_n_b;
% end
% 
% r_mag = R_surface(Earth_R_long, Earth_R_short, P1(1,:));%*(1+tan(wr)^2);


%% 역 연산
%R_surface(Earth_R_long, Earth_R_short, P1(1,:));%sqrt((1+tan(P1(1,:)).^2)./(Earth_R_short^2 + (Earth_R_long^2)*tan(P1(1,:)).^2))*Earth_R_short*Earth_R_long + P1(3,i);
%R_vector = r_mag.*[cos(P1(1,:)).*cos(P1(2,:)); cos(P1(1,:)).*sin(P1(2,:)); sin(P1(1,:))]; %P1(3,:).*
surf_r = R_surface(Earth_R_long, Earth_R_short, P1(1,:)) + P1(3,:);
r = [surf_r.*cos(P1(1,:)).*cos(P1(2,:)); surf_r.*cos(P1(1,:)).*sin(P1(2,:)); surf_r.*sin(P1(1,:))];%g = GM*1./r_mag.^2;
g0 = 9.780318 * (1 + 5.3024e-3 .* sin(P1(1,:)).*sin(P1(1,:)) - 5.9e-6 .* sin(2*P1(1,:)) .* sin(2*P1(1,:)));
g = g0./(1+P1(3,:)/Earth_R_long);
w_n_ie = [Earth_Omega*cos(P1(1,:)); zeros(size(P1(1,:))); -Earth_Omega*sin(P1(1,:))];
g_n_l = [zeros(size(g)); zeros(size(g)); g] - cross(w_n_ie, cross(w_n_ie, r)); %애초에 중력을 일케 쓰면 무조건 로컬

delta_L=zeros(1,size(P1,2));delta_l=zeros(1,size(P1,2));delta_h=zeros(1,size(P1,2));RE=zeros(1,size(P1,2));RN=zeros(1,size(P1,2));
for i=1:size(P1,2)-1
    RE(i) = Earth_R_long/sqrt(1-(e*sin(P1(1,i)))^2);
    RN(i) = Earth_R_long*(1-e^2)/(1-(e*sin(P1(1,i)))^2)^(3/2);
    delta_L(i) = (P1(1,i+1) - P1(1,i))/dt; 
    delta_l(i) = (P1(2,i+1) - P1(2,i))/dt; 
    delta_h(i) = (P1(3,i+1) - P1(3,i))/dt; 
end


w_n_en = [delta_l.*cos(P1(1,:)); -delta_L; -delta_l.*sin(P1(1,:))];
w_n_in = w_n_ie + w_n_en;
w_b_in = zeros(3, size(P1,2));
for i=1:size(P1,2)
    w_b_in(:,i) = qC(:,:,i)' * w_n_in(:,i);
end

% delta_qC = zeros(3,3,size(P1,2));
% Omega_b_nb = zeros(3,3,size(P1,2));
% w_b_nb = zeros(3, size(P1,2));
% for i=1:size(P1,2)-1
%     delta_qC(:,:,i) = (qC(:,:,i+1) - qC(:,:,i))/dt;
% end
% for i=1:size(P1,2)
%     Omega_b_nb(:,:,i) = qC(:,:,i)'*delta_qC(:,:,i);
%     w_b_nb(:,i) = w_nb;%[Omega_b_nb(3,2,i); Omega_b_nb(1,3,i); Omega_b_nb(2,1,i)];   
% end
w_b_ib = w_b_in + w_b_nb; %자이로 센서 데이터

V_N = delta_L.*(P1(3,:)+RN);
V_E = delta_l.*(P1(3,:)+RE).*cos(P1(1,:));
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
   f_b(:,i) = qC(:,:,i)'*f_n(:,i);
end
% inv_f_n = f_n;
% inv_qC = qC;
% inv_q = q;
% inv_w_b_nb = w_b_nb;
% inv_w_b_in = w_b_in;
% inv_w_n_en = w_n_en;
% inv_w_n_ie = w_n_ie;
 inv_V_n = V_n;
 inv_P = P1;
% inv_dL = delta_L;
% inv_dl = delta_l;
% inv_dh = delta_h;

% [gBias, gARW, aBias, aVRW] = IMUspec('HG1700AG58');
% 
% gRandomwalk = zeros(3,size(w_b_ib,2));
% aRandomwalk = zeros(3,size(w_b_ib,2));
% for i=1:size(w_b_ib,2)
%    gRandomwalk(:,i) = gARW/60*sqrt(i*dt)*(rand(size(w_b_ib(:,i)))-0.5)*2; %(-1, 1) randomwalk?
%    aRandomwalk(:,i) = aVRW;
% end
% 
% %error term
% w_b_ib = w_b_ib + (gBias/3600)*randn() + gRandomwalk;
% f_b = f_b  + aBias*g/1000*randn() + aRandomwalk;


save Qbody_data.mat f_b w_b_ib init_P init_q init_V_n inv_V_n inv_P; %inv_dL inv_dl inv_dh inv_f_n inv_qC inv_q inv_w_b_nb inv_w_b_in inv_w_n_en inv_w_n_ie 

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




