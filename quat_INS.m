quat_inverseINS;
clc;
clear;

load Qbody_data.mat; % load f_b, w_b_ib, init_P, init_q, init_V_n

dt = 0.01; %100Hz
Earth_Omega = 7.292115e-5;
Earth_R_short = 6356752.3142;
Earth_R_long = 6378137.0;
GM = 398600.4418*1000^3;
[ex, ey, ez] = ellipsoid(0,0,0,Earth_R_long,Earth_R_long,Earth_R_short,20);

%% 역 연산

% 초기값
P = init_P'; % 위도(L), 경도(l), 중심으로부터 거리(r)

q = init_q;
qC = quat2DCM(q(:,1));
V_n = init_V_n;

% INS

for i=1:size(f_b,2)
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
    wb = [0, 10, 0];
    q(:,i+1) = quat_update(q(:,i),w_b_nb(:,i),dt);
    q(:,i+1) = q(:,i+1)/norm(q(:,i+1)); %normalization
    qC(:,:,i+1) = quat2DCM(q(:,i+1)); %DCM
    
    f_n(:,i) = qC(:,:,i) * f_b(:,i);
    
    %r_mag(i) = R_surface(Earth_R_long, Earth_R_short, P(1,i));%sqrt((1+tan(P(1,i))^2)./(Earth_R_short^2 + (Earth_R_long^2)*tan(P(1,i))^2))*Earth_R_short*Earth_R_long;
    %r = r_mag(i)*[cos(P(1,i))*cos(P(2,i)); cos(P(1,i))*sin(P(2,i)); sin(P(1,i))];
    r = [P(3,i)*cos(P(1,i))*cos(P(2,i)), P(3,i)*cos(P(1,i))*sin(P(2,i)), P(3,i)*sin(P(1,i))]';
    g0 = 9.780318 * (1 + 5.3024e-3 .* sin(P(1,i)).*sin(P(1,i)) - 5.9e-6 .* sin(2*P(1,i)) .* sin(2*P(1,i)));
    g = g0/(1+(P(3,i)-Earth_R_long)/Earth_R_long);
    g_n_l(:,i) = [zeros(size(g)); zeros(size(g)); g] - cross(w_n_ie(:,i), cross(w_n_ie(:,i), r));
    
    delta_V_n = f_n(:,i) - cross(2*w_n_ie(:,i) + w_n_en(:,i), V_n(:,i)) + g_n_l(:,i);
    V_n(:,i+1) = V_n(:,i) + delta_V_n*dt;
    P(:,i+1) = P(:,i) + [delta_L; delta_l; delta_h]*dt;
    sL(i) = delta_L; sl(i) = delta_l; sh(i) = delta_h;
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
