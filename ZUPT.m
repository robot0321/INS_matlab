clc;
clear all;

load body_data.mat; % load f_b, w_b_ib, init_P, init_C_n_b, init_V_n, true_P, true_V, true_C_n_b

dt = 0.01; %100Hz
f_b_ib = f_b;
N = size(f_b,2)-1;

Earth_Omega = 7.292115e-5;
Earth_R_long = 6378137.0;
e = 0.0818191908425;

P = init_P-[0;0;Earth_R_long];
C_n_b = init_C_n_b;
V_n = init_V_n;
x = zeros(15,1);
PP = 100*eye(15,15);

% Kalman_real
% x, PP, Q, R, z, C_n_b, f_b_ib, w_n_ie, w_n_en, v_n_eb, Lb, hb, dt, l_b_ba
    
for i=1:N
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
%     delta_C_n_b(:,:,i) = C_n_b(:,:,i)*Omega_b_nb(:,:,i);
%     Omega_b_ib = angularV2M(w_b_ib(:,i));
%     Omega_n_in = angularV2M(w_n_in);
%     delta_C_n_b(:,:,i) = C_n_b(:,:,i)*Omega_b_ib - Omega_n_in*C_n_b(:,:,i);
    Omega_b_ib = angularV2M(w_b_ib(:,i));
    Omega_n_ie = angularV2M(w_n_ie(:,i));
    Omega_n_en = angularV2M(w_n_en(:,i));
    delta_C_n_b(:,:,i) = C_n_b(:,:,i)*Omega_b_ib - (Omega_n_ie+Omega_n_en)*C_n_b(:,:,i);
    C_n_b(:,:,i+1) = C_n_b(:,:,i) + delta_C_n_b(:,:,i)*dt;
    

    f_n(:,i) = (C_n_b(:,:,i) + C_n_b(:,:,i+1))/2 * f_b(:,i);
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
    
for i=1:1
    %%
    Lb = P(i,1);
    hb = P(i,3)-Earth_R_long;
  
    RE = Earth_R_long/sqrt(1-(e*sin(Lb))^2);
    RN = Earth_R_long*(1-e^2)/(1-(e*sin(Lb))^2)^(3/2);
    g0 = 9.7803253359*(1+0.001931853*sin(Lb)^2)/sqrt(1-(e*sin(Lb))^2);
    r_e_eS = RE*sqrt(1-e*(2-e)*sin(Lb)^2);
    
    F11 = -angularV2M(w_n_ie(:,i) + w_n_en(:,i));
    F12 = [0, -1/(RE + hb), 0;
           1/(RN + hb), 0, 0;
           0, tan(Lb)/(RE+hb), 0];
    F13 = [Earth_Omega*sin(Lb), 0, v_n_eb(2)/(RE+hb)^2;
           0, 0, -v_n_eb(1)/(RN+hb)^2;
           Earth_Omega*cos(Lb) +v_n_eb(2)/(RE+hb)/(cos(Lb))^2, 0, -v_n_eb(2)*tan(Lb)/(RE+hb)^2];
    F21 = -angularV2M(C_n_b*f_b_ib);
    F22 = [v_n_eb(3)/(RN+hb), -2*v_n_eb(2)*tan(Lb)/(RE+hb) -2*Earth_Omega*sin(Lb), v_n_eb(1)/(RN+hb);
           v_n_eb(2)*tan(Lb)/(RE+hb) +2*Earth_Omega*sin(Lb), (v_n_eb(1)*tan(Lb) +v_n_eb(3))/(RE+hb), v_n_eb(2)/(RE+hb) +2*Earth_Omega*cos(Lb);
           -2*v_n_eb(1)/(RN+hb), -2*v_n_eb(2)/(RE+hb) -2*Earth_Omega*cos(Lb), 0];
    F23 = [-(v_n_eb(2)^2)*(sec(Lb)^2)/(RE+hb) -2*v_n_eb(2)*Earth_Omega*cos(Lb), 0, v_n_eb(2)^2*tan(Lb)/(RE+hb)^2 -v_n_eb(1)*v_n_eb(3)/(RN+hb)^2;
           v_n_eb(1)*v_n_eb(2)*sec(Lb)^2/(RE+hb) +2*v_n_eb(1)*Earth_Omega*cos(Lb) -2*v_n_eb(3)*Earth_Omega*sin(Lb), 0, -(v_n_eb(1)*v_n_eb(2)*tan(Lb) +v_n_eb(2)*v_n_eb(3))/(RE+hb)^2;
           2*v_n_eb(2)*Earth_Omega*sin(Lb), 0, v_n_eb(2)^2/(RE+hb)^2 +v_n_eb(1)^2/(RN+hb)^2 -2*g0/r_e_eS];
    F32 = [1/(RN+hb), 0, 0;
           0, 1/((RE+hb)*cos(Lb)), 0;
           0, 0, -1];
    F33 = [0, 0, -v_n_eb(1)/(RN+hb)^2;
           v_n_eb(2)*sin(Lb)/((RE+hb)*cos(Lb)^2), 0, v_n_eb(2)/(RE+hb)^2/cos(Lb);
           0, 0, 0];
    
    zr3 = zeros(3,3);
    F = [F11, F12, F13, zr3, C_n_b;
         F21, F22, F23, C_n_b, zr3;
         zr3, F32, F33, zr3, zr3;
         zr3, zr3, zr3, zr3, zr3;
         zr3, zr3, zr3, zr3, zr3];
    A = eye(15) + F*dt;

    Q = zeros(15,15);
    Q(1:3,1:3) = eye(3)*0.01;
    Q(4:6,4:6) = eye(3)*3e-3;
    Q(10:12,10:12) = eye(3)*0.01;
    Q(13:15,13:15) = eye(3)*0.001;
    
    R = 100*eye(15,15)*i*dt;
    
    H = [zr3, -eye(3), zr3, zr3, zr3]; % 
    z = -[0;0;0];
    
    xp = A*x;
    Pp = A*PP*A' + Q;
    K = Pp*H'/(H*Pp*H' + R);
    x_k1 = xp + K*(z - H*xp);
    P_k1 = Pp - K*H*Pp;
    
    PP = Pp;
    x = xp;
    
    
end
