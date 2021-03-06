clc;
clear all;
close all;

load body_data.mat; % load f_b, w_b_ib, init_P, init_C_n_b, init_V_n, true_P, true_V, true_C_n_b

dt = 0.01; %100Hz
f_b_ib = f_b;
N = size(f_b,2)-1;

Earth_Omega = 7.292115e-5;
Earth_R_long = 6378137.0;
e = 0.0818191908425;

P = init_P'-[0;0;Earth_R_long];
C_n_b = init_C_n_b;
R_vector = [(P(3,1)+Earth_R_long)*cos(P(1,1))*cos(P(2,1)); (P(3,1)+Earth_R_long)*cos(P(1,1))*sin(P(2,1)); (P(3,1)+Earth_R_long)*sin(P(1,1))];
w_n_ie(:,1) = [Earth_Omega*cos(P(1,1)); 0; -Earth_Omega*sin(P(1,1))];
v_n_eb(:,1) = init_V_n;% - Earth_Omega*norm(R_vector)*cos(P(1,1));
    
x_k = 0.1*ones(15,1);
P_k = 1000*eye(15,15);

% Kalman_real
% x, PP, Q, R, z, C_n_b, f_b_ib, w_n_ie, w_n_en, v_n_eb, Lb, hb, dt, l_b_ba
    
for i=1:N
    Lb = P(1,i);
    hb = P(3,i);
  
    RE(i) = Earth_R_long/sqrt(1-(e*sin(Lb))^2);
    RN(i) = Earth_R_long*(1-e^2)/(1-(e*sin(Lb))^2)^(3/2);
    
    % N E D
%     delta_L = V_n(1,i)/(P(3,i));
%     delta_l = V_n(2,i)/(P(3,i))/cos(P(1,i));
%     delta_h = -V_n(3,i);
    
    %w_n_en(:,i) = [delta_l*cos(P(1,i)); -delta_L; -delta_l*sin(P(1,i))];
    w_n_en(:,i) = [v_n_eb(2,i)/(RE(i)+hb); -v_n_eb(1,i)/(RN(i)+hb); -v_n_eb(2,i)*tan(P(1,i))/(RE(i)+hb)];
    w_n_ie(:,i) = [Earth_Omega*cos(P(1,i)); 0; -Earth_Omega*sin(P(1,i))];
    w_n_in(:,i) = w_n_ie(:,i) + w_n_en(:,i);
    w_b_in(:,i) = C_n_b(:,:,i)' * w_n_in(:,i);
    w_b_nb(:,i) = w_b_ib(:,i) - w_b_in(:,i);
    
    %Attitude
    Omega_b_nb(:,:,i) = angularV2M(w_b_nb(:,i));%[0, -w_b_nb(3), w_b_nb(2); w_b_nb(3), 0, -w_b_nb(1); -w_b_nb(2), w_b_nb(1), 0];
    delta_C_n_b = C_n_b(:,:,i)*Omega_b_nb(:,:,i);
%     Omega_b_ib = angularV2M(w_b_ib(:,i));
%     Omega_n_in = angularV2M(w_n_in);
%     delta_C_n_b(:,:,i) = C_n_b(:,:,i)*Omega_b_ib - Omega_n_in*C_n_b(:,:,i);
    Omega_b_ib = angularV2M(w_b_ib(:,i));
    Omega_n_ie = angularV2M(w_n_ie(:,i));
    Omega_n_en = angularV2M(w_n_en(:,i));
    %delta_C_n_b = C_n_b(:,:,i)*Omega_b_ib - (Omega_n_ie+Omega_n_en)*C_n_b(:,:,i);
    C_n_b(:,:,i+1) = C_n_b(:,:,i) + delta_C_n_b*dt;
    
    %specific force frame transformation
    f_n(:,i) = (C_n_b(:,:,i) + C_n_b(:,:,i+1))/2 * f_b(:,i);
    
    
    %r_mag(i) = R_surface(Earth_R_long, Earth_R_short, P(1,i));%sqrt((1+tan(P(1,i))^2)./(Earth_R_short^2 + (Earth_R_long^2)*tan(P(1,i))^2))*Earth_R_short*Earth_R_long;
    %r = r_mag(i)*[cos(P(1,i))*cos(P(2,i)); cos(P(1,i))*sin(P(2,i)); sin(P(1,i))];
    r = [P(3,i)*cos(P(1,i))*cos(P(2,i)), P(3,i)*cos(P(1,i))*sin(P(2,i)), P(3,i)*sin(P(1,i))]';
    %g0 = 9.780318 * (1 + 5.3024e-3 .* sin(P(1,i)).*sin(P(1,i)) - 5.9e-6 .* sin(2*P(1,i)) .* sin(2*P(1,i)));
    g0 = 9.7803253359*(1+0.001931853*sin(Lb)^2)/sqrt(1-(e*sin(Lb))^2);
    r_e_eS = RE(i)*sqrt(1-e*(2-e)*sin(Lb)^2);
    g = g0/(1+P(3,i)/Earth_R_long);
    g_n_l(:,i) = [zeros(size(g)); zeros(size(g)); g] - cross(w_n_ie(:,i), cross(w_n_ie(:,i), r));
    
    %velocity update
    delta_v_n_eb = f_n(:,i) - cross(2*w_n_ie(:,i) + w_n_en(:,i), v_n_eb(:,i)) + g_n_l(:,i);
    v_n_eb(:,i+1) = v_n_eb(:,i) + delta_v_n_eb*dt;
    
    P(3,i+1) = P(3,i) - (v_n_eb(3,i+1)+v_n_eb(3,i))/2*dt; %h
    P(1,i+1) = P(1,i) + (v_n_eb(1,i+1)/(RN(i)+P(3,i+1)) + v_n_eb(1,i)/(RN(i)+hb))/2*dt; %L
    RE(i+1) = Earth_R_long/sqrt(1-(e*sin(P(1,i+1)))^2);
    P(2,i+1) = P(2,i) + (v_n_eb(2,i+1)/(RE(i+1)+P(3,i+1))/cos(P(1,i+1)) + v_n_eb(2,i)/(RE(i)+hb)/cos(Lb))/2*dt; %lambda
    
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
    Q(1:3,1:3) = eye(3)*0.01;
    Q(4:6,4:6) = eye(3)*3e-3;
    Q(10:12,10:12) = eye(3)*0.01;
    Q(13:15,13:15) = eye(3)*0.001;
    
    R = 100*eye(3,3); %선정 어뜨케 하지
    
    H = [zr3, -eye(3), zr3, zr3, zr3]; % 
    z = -[0;0;0];
    
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
    
    if i==1000
        pause(1);
    end
    
end

figure(1);
subplot(3,1,1);
hold on;
plot(P(1,1:end-1)*180/pi)
plot(true_P(1,1:end-1)*180/pi)
subplot(3,1,2);
hold on;
plot(P(2,1:end-1)*180/pi)
plot(true_P(2,1:end-1)*180/pi)
subplot(3,1,3);
hold on;
plot(P(3,1:end-1))
plot(true_P(3,1:end-1)-Earth_R_long)

figure(2);
subplot(3,1,1);
hold on;
plot(v_n_eb(1,1:end-1)*180/pi)
plot(true_V(1,1:end-1)*180/pi)
subplot(3,1,2);
hold on;
plot(v_n_eb(2,1:end-1)*180/pi)
plot(true_V(2,1:end-1)*180/pi)
subplot(3,1,3);
hold on;
plot(v_n_eb(3,1:end-1))
plot(true_V(3,1:end-1))