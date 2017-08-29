quat_inverseINS;
clc;
clear;

load Qbody_data.mat; % load f_b, w_b_ib, init_P, init_q, init_V_n


% % add noise
% 
% % specification: MPU9250
% f_b = f_b + 0.08*randn(size(f_b)) + 0.01; % noise + bias
% w_b_ib = w_b_ib + 0.1*pi/180*randn(size(w_b_ib)) + 0.001; % noise + bias

% specification: HG1700
% f_b = f_b + 0.001*randn(size(f_b));  % noise + bias
% w_b_ib = w_b_ib + 0.001*randn(size(w_b_ib)); % noise + bias

% %add noise w_b



dt = 0.01; %100Hz
Earth_Omega = 7.292115e-5;
Earth_R_short = 6356752.3142;
Earth_R_long = 6378137.0;
GM = 398600.4418*1000^3;
[ex, ey, ez] = ellipsoid(0,0,0,Earth_R_long,Earth_R_long,Earth_R_short,20);
e = 0.0818191908425;
    
x_k = 0.1*ones(15,1);
P_k = 100*eye(15,15);


%% 역 연산

% 초기값
P = init_P'; % 위도(L), 경도(l), 중심으로부터 거리(r)

q = init_q;
qC = quat2DCM(q(:,1));
V_n = init_V_n;

C_n_b = quat2DCM(init_q);
R_vector = [(P(3,1)+Earth_R_long)*cos(P(1,1))*cos(P(2,1)); (P(3,1)+Earth_R_long)*cos(P(1,1))*sin(P(2,1)); (P(3,1)+Earth_R_long)*sin(P(1,1))];
w_n_ie(:,1) = [Earth_Omega*cos(P(1,1)); 0; -Earth_Omega*sin(P(1,1))];
v_n_eb(:,1) = init_V_n;% - Earth_Omega*norm(R_vector)*cos(P(1,1));

% INS

for i=1:size(f_b,2)-1
    RE(i) = Earth_R_long/sqrt(1-(e*sin(P(1,i)))^2);
    
    delta_L = V_n(1,i)/(P(3,i));
    delta_l = V_n(2,i)/(P(3,i))/cos(P(1,i));
    delta_h = -V_n(3,i);

    w_n_en(:,i) = [delta_l*cos(P(1,i)); -delta_L; -delta_l*sin(P(1,i))];
    w_n_ie(:,i) = [Earth_Omega*cos(P(1,i)); 0; -Earth_Omega*sin(P(1,i))];
    w_n_in(:,i) = w_n_ie(:,i) + w_n_en(:,i);
    w_b_in(:,i) = qC(:,:,i)' * w_n_in(:,i);
    w_b_nb(:,i) = w_b_ib(:,i) - w_b_in(:,i);
    
    
    %Omega_b_nb(:,:,i) = angularV2M(w_b_nb(:,i));%[0, -w_b_nb(3), w_b_nb(2); w_b_nb(3), 0, -w_b_nb(1); -w_b_nb(2), w_b_nb(1), 0];
    %delta_qC(:,:,i) = qC(:,:,i)*Omega_b_nb(:,:,i);
%     Omega_b_ib = angularV2M(w_b_ib(:,i));
%     Omega_n_in = angularV2M(w_n_in);
%     delta_qC(:,:,i) = qC(:,:,i)*Omega_b_ib - Omega_n_in*qC(:,:,i);
    %qC(:,:,i+1) = qC(:,:,i) + delta_qC(:,:,i)*dt;
    %qC(:,:,i+1) = DCM_ortho_normal_compensation(qC(:,:,i+1));
    q(:,i+1) = quat_update(q(:,i),w_b_nb(:,i),dt);
    q(:,i+1) = q(:,i+1)/norm(q(:,i+1)); %normalization
    qC(:,:,i+1) = quat2DCM(q(:,i+1)); %DCM
    
    f_n(:,i) = qC(:,:,i) * f_b(:,i);
    
    %r_mag(i) = R_surface(Earth_R_long, Earth_R_short, P(1,i));%sqrt((1+tan(P(1,i))^2)./(Earth_R_short^2 + (Earth_R_long^2)*tan(P(1,i))^2))*Earth_R_short*Earth_R_long;
    %r = r_mag(i)*[cos(P(1,i))*cos(P(2,i)); cos(P(1,i))*sin(P(2,i)); sin(P(1,i))];
    r = [P(3,i)*cos(P(1,i))*cos(P(2,i)), P(3,i)*cos(P(1,i))*sin(P(2,i)), P(3,i)*sin(P(1,i))]';
    %g0 = 9.7803253359*(1+0.001931853*sin(P(1,i))^2)/sqrt(1-(e*sin(P(1,i)))^2);
    g0 = 9.780318 * (1 + 5.3024e-3 .* sin(P(1,i)).*sin(P(1,i)) - 5.9e-6 .* sin(2*P(1,i)) .* sin(2*P(1,i)));
    r_e_eS = RE(i)*sqrt(1-e*(2-e)*sin(P(1,i))^2);
    g = g0/(1+(P(3,i)-Earth_R_long)/Earth_R_long);
    g_n_l(:,i) = [zeros(size(g)); zeros(size(g)); g] - cross(w_n_ie(:,i), cross(w_n_ie(:,i), r));
    

    delta_V_n = f_n(:,i) - cross(2*w_n_ie(:,i) + w_n_en(:,i), V_n(:,i)) + g_n_l(:,i);
    V_n(:,i+1) = V_n(:,i) + delta_V_n*dt;
    P(:,i+1) = P(:,i) + [delta_L; delta_l; delta_h]*dt;
    sL(i) = delta_L; sl(i) = delta_l; sh(i) = delta_h;
    
    
    
    %추가
    Lb = P(1,i);
    hb = P(3,i);
    RE(i) = Earth_R_long/sqrt(1-(e*sin(Lb))^2);
    RN(i) = Earth_R_long*(1-e^2)/(1-(e*sin(Lb))^2)^(3/2);
    v_n_eb(:,i) = V_n(:,i);
    v_n_eb(:,i+1) = V_n(:,i+1);
    C_n_b(:,:,i+1) = qC(:,:,i+1);
    
    % Kalman
    % x = (d_att3, d_vel3, d_pos3, ba3, bg3)
    F11 = -angularV2M(w_n_ie(:,i) + w_n_en(:,i));
    F12 = [0, -1/(RE(i) + hb), 0;
           1/(RN(i) + hb), 0, 0;
           0, tan(Lb)/(RE(i)+hb), 0];
    F13 = [Earth_Omega*sin(Lb), 0, v_n_eb(2,i)/(RE(i)+hb)^2;
           0, 0, -v_n_eb(1,i)/(RN(i)+hb)^2;
           Earth_Omega*cos(Lb)+v_n_eb(2,i)/(RE(i)+hb)/(cos(Lb))^2, 0, -v_n_eb(2,i)*tan(Lb)/(RE(i)+hb)^2];
    F21 = -angularV2M(f_n(:,i));
    F22 = [v_n_eb(3,i)/(RN(i)+hb), -2*v_n_eb(2,i)*tan(Lb)/(RE(i)+hb)-2*Earth_Omega*sin(Lb), v_n_eb(1,i)/(RN(i)+hb);
           v_n_eb(2,i)*tan(Lb)/(RE(i)+hb)+2*Earth_Omega*sin(Lb), (v_n_eb(1,i)*tan(Lb)+v_n_eb(3,i))/(RE(i)+hb), v_n_eb(2,i)/(RE(i)+hb)+2*Earth_Omega*cos(Lb);
           -2*v_n_eb(1,i)/(RN(i)+hb), -2*v_n_eb(2,i)/(RE(i)+hb)-2*Earth_Omega*cos(Lb), 0];
    F23 = [-(v_n_eb(2,i)^2)*(sec(Lb)^2)/(RE(i)+hb)-2*v_n_eb(2,i)*Earth_Omega*cos(Lb), 0, v_n_eb(2,i)^2*tan(Lb)/(RE(i)+hb)^2-v_n_eb(1,i)*v_n_eb(3,i)/(RN(i)+hb)^2;
           v_n_eb(1,i)*v_n_eb(2,i)*sec(Lb)^2/(RE(i)+hb)+2*v_n_eb(1,i)*Earth_Omega*cos(Lb)-2*v_n_eb(3,i)*Earth_Omega*sin(Lb), 0, -(v_n_eb(1,i)*v_n_eb(2,i)*tan(Lb)+v_n_eb(2,i)*v_n_eb(3,i))/(RE(i)+hb)^2;
           2*v_n_eb(2,i)*Earth_Omega*sin(Lb), 0, v_n_eb(2,i)^2/(RE(i)+hb)^2+v_n_eb(1,i)^2/(RN(i)+hb)^2-2*g0/r_e_eS];
    F32 = [1/(RN(i)+hb), 0, 0;
           0, 1/((RE(i)+hb)*cos(Lb)), 0;
           0, 0, -1];
    F33 = [0, 0, -v_n_eb(1,i)/(RN(i)+hb)^2;
           v_n_eb(2,i)*sin(Lb)/((RE(i)+hb)*cos(Lb)^2), 0, v_n_eb(2,i)/(RE(i)+hb)^2/cos(Lb);
           0, 0, 0];
    
    zr3 = zeros(3,3);
    F = [F11, F12, F13, zr3, C_n_b(:,:,i+1);
         F21, F22, F23, C_n_b(:,:,i+1), zr3;
         zr3, F32, F33, zr3, zr3;
         zr3, zr3, zr3, zr3, zr3;
         zr3, zr3, zr3, zr3, zr3];
    A = eye(15) + F*dt;

    Q = zeros(15,15);
    Q(1:3,1:3) = eye(3)*1e-2;
    Q(4:6,4:6) = eye(3)*3e-3;
    Q(10:12,10:12) = eye(3)*1e-2;
    Q(13:15,13:15) = eye(3)*1e-3;
    
    R = 100*eye(3);%diag([3/Earth_R_long, 3/Earth_R_long, 20]); %선정 어뜨케 하지 크게(10000000000*eye(3)) 잡으면 잘된다?
    
%     %ZUPT
%     H = [zr3, -eye(3), zr3, zr3, zr3]; % 
%     z = -(v_n_eb(:,i+1) - [0;0;0]);
    
%     %CUPT
%     H = [zr3, zr3, -eye(3), zr3, zr3]; % 
%     z = -(P(:,i+1) - init_P');    
    
    %GPS
    GPS_sd = [3/Earth_R_long, 3/Earth_R_long, 20]'/5; %3m, 3m, 20m 오차
    GPS_P = inv_P(:,i+1);% + GPS_sd.*randn(3,1);
    H = [zr3, zr3, -eye(3), zr3, zr3]; % 
    z = -(P(:,i+1) - GPS_P);
    GG(:,i) = GPS_P;
    
    
    xp = A*x_k;
    Pp = A*P_k*A' + Q;
    K = Pp*H'/(H*Pp*H' + R);
    x_k1 = xp + K*(z - H*xp);
    P_k1 = Pp - K*H*Pp;
    
    P_k = P_k1;
    x_k = x_k1;
    
    C_n_b(:,:,i+1) = (eye(3) - angularV2M(x_k(1:3)))*C_n_b(:,:,i+1);
    v_n_eb(:,i+1) = v_n_eb(:,i+1) - x_k(4:6);
    P(:,i+1) = P(:,i+1) - x_k(7:9);
    
    qC(:,:,i+1) = C_n_b(:,:,i+1);
    %추가
    
    V_n(:,i+1) = v_n_eb(:,i+1);
    
    sx_k(:,i) = x_k;
end

figure(1);
subplot(3,1,1);
hold on;
plot(P(1,1:end-1)*180/pi)
subplot(3,1,2);
hold on;
plot(P(2,1:end-1)*180/pi)
subplot(3,1,3);
hold on;
plot(P(3,1:end-1))

figure(2);
%surf(ex, ey, ez);
%axis equal
hold on;
plot3(P(3,1:end-1).*cos(P(1,1:end-1)).*cos(P(2,1:end-1)), P(3,1:end-1).*cos(P(1,1:end-1)).*sin(P(2,1:end-1)), P(3,1:end-1).*sin(P(1,1:end-1)), 'LineWidth', 2, 'Marker', '.');
figure(3);
hold on;
plot(P(1,1:end-1)*180/pi,P(2,1:end-1)*180/pi);
% 
% figure();
% attitude_visualize(qC,qC);

% figure(3);
% subplot(3,1,1);
% hold on;
% plot(V_n(1,1:end-1))
% subplot(3,1,2);
% hold on;
% plot(V_n(2,1:end-1))
% subplot(3,1,3);
% hold on;
% plot(V_n(3,1:end-1))
