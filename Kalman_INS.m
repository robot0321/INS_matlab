%inverseINS;
close all
clc;
clear all;

load body_data.mat; % load f_b, w_b_ib, init_P, init_C_n_b, init_V_n, true_P, true_V, true_C_n_b

dt = 0.01; %100Hz
Earth_Omega = 7.292115e-5;
Earth_R_short = 6356752.3142;
Earth_R_long = 6378137.0;
GM = 398600.4418*1000^3;
[ex, ey, ez] = ellipsoid(0,0,0,Earth_R_long,Earth_R_long,Earth_R_short,20);

%% 역 연산

% 초기값
% NED NED NED NED NED
P = init_P'; % 위도(L), 경도(l), 중심으로부터 거리(r)
C_n_b = init_C_n_b;
V_n = init_V_n;


err_x = zeros(15,1); %state vector means 'error'
Pv = 100*eye(15,15);
Q = 0.01*eye(15,15);
% Q(1:3,1:3) = eye(3);
% Q(4:6,4:6) = eye(3);
Q(7:9,7:9) = zeros(3);
% Q(10:12,10:12) = eye(3);
% Q(13:15,13:15) = eye(3);
R = 10*eye(3);

% INS

for i=1:size(f_b,2)-1
    if i==85
        pause(3);
    end
    % N E D
    delta_L = V_n(1,i)/(P(3,i));
    delta_l = V_n(2,i)/(P(3,i))/cos(P(1,i));
    delta_h = -V_n(3,i);

    w_n_en(:,i) = [delta_l*cos(P(1,i)); -delta_L; -delta_l*sin(P(1,i))];
    w_n_ie(:,i) = [Earth_Omega*cos(P(1,i)); 0; -Earth_Omega*sin(P(1,i))];
    w_n_in(:,i) = w_n_ie(:,i) + w_n_en(:,i);
    w_b_in(:,i) = C_n_b(:,:,i)' * w_n_in(:,i);
    w_b_nb(:,i) = w_b_ib(:,i) - w_b_in(:,i);
    
    R_vector = [P(3,i).*cos(P(1,i)).*cos(P(2,i)); P(3,i).*cos(P(1,i)).*sin(P(2,i)); P(3,i).*sin(P(1,i))];
    v_n_ie = Earth_Omega*norm(R_vector)*cos(P(1,i));
    v_n_eb = V_n(:,i) - v_n_ie;
    
    Omega_b_nb(:,:,i) = angularV2M(w_b_nb(:,i));%[0, -w_b_nb(3), w_b_nb(2); w_b_nb(3), 0, -w_b_nb(1); -w_b_nb(2), w_b_nb(1), 0];
    delta_C_n_b(:,:,i) = C_n_b(:,:,i)*Omega_b_nb(:,:,i);
%     Omega_b_ib = angularV2M(w_b_ib(:,i));
%     Omega_n_in = angularV2M(w_n_in);
%     delta_C_n_b(:,:,i) = C_n_b(:,:,i)*Omega_b_ib - Omega_n_in*C_n_b(:,:,i);
    C_n_b(:,:,i+1) = C_n_b(:,:,i) + delta_C_n_b(:,:,i)*dt;
    

    f_n(:,i) = C_n_b(:,:,i) * f_b(:,i);
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
    
%     %************** Measurement
%     noiseP = zeros(3,1);%3*randn(3,1);
%     noiseV = 0.004*randn(3,1);
% %     true_z = true_V(:,i+1); 
% %     propa_z = V_n(:,i+1); 
%     true_z = true_P(:,i+1)-[0,0,Earth_R_long]'+noiseP;
%     propa_z = P(:,i+1)-[0,0,Earth_R_long]'+noiseP; 
%     
%     err_z = propa_z - true_z;
%     %************** Kalman    
%     [err_x(:,i+1), Pv] = Kal_AttVelPos(err_x(:,i), Pv, Q, R, err_z, C_n_b(:,:,i), f_b(:,i), w_n_ie(:,i), w_n_en(:,i), v_n_eb, P(1,i), P(3,i)-Earth_R_long, dt, [0,0,0]');
%     
%     %************** Kalman solution corretion
%     C_n_b(:,:,i+1) = (eye(3) - angularV2M(err_x(1:3,i+1)))*C_n_b(:,:,i+1);
%     %C_n_b(:,:,i+1) = DCM_ortho_normal_compensation(C_n_b(:,:,i+1));
%     V_n(:,i+1) = V_n(:,i+1) - err_x(4:6,i+1); % 원래는 v_n_eb = v_n_eb - x(4:6)';이고 V_n계산에서도 v_n_ie가 P가 바뀜에 따라 바뀌어서 바꿔줘야 되지만 근사.
% %     P(:,i+1) = P(:,i+1) - err_x(7:9,i+1);
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

