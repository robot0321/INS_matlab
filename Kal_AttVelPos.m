function [x_k1, P_k1] = Kal_AttVelPos(x, PP, Q, R, z, C_n_b, f_b_ib, w_n_ie, w_n_en, v_n_eb, Lb, hb, dt, l_b_ba)
    % x = [att3; vel3; pos3; acc_bias; gyro_bias];
    % NED frame means, v(1) = v_N, v(2) = v_E, v(3) = v_D;
    persistent RE RN r_e_eS F11 F12 F13 
    
    Earth_Omega = 7.292115e-5;
    Earth_R_long = 6378137.0;
    e = 0.0818191908425;
    
    RE = Earth_R_long/sqrt(1-(e*sin(Lb))^2);
    RN = Earth_R_long*(1-e^2)/(1-(e*sin(Lb))^2)^(3/2);
    g0 = 9.7803253359*(1+0.001931853*sin(Lb)^2)/sqrt(1-(e*sin(Lb))^2);
    r_e_eS = RE*sqrt(1-e*(2-e)*sin(Lb)^2);
    
    F11 = -angularV2M(w_n_ie + w_n_en);
    F12 = [0, -1/(RE + hb), 0;
           1/(RN + hb), 0, 0;
           0, tan(Lb)/(RE+hb), 0];
    F13 = [Earth_Omega*sin(Lb), 0, v_n_eb(2)/(RE+hb)^2;
           0, 0, -v_n_eb(1)/(RN+hb)^2;
           Earth_Omega*cos(Lb)+v_n_eb(2)/(RE+hb)/(cos(Lb))^2,0,-v_n_eb(2)*tan(Lb)/(RE+hb)^2];
    F21 = -angularV2M(C_n_b*f_b_ib);
    F22 = [v_n_eb(3)/(RN+hb), -2*v_n_eb(2)*tan(Lb)/(RE+hb)-2*Earth_Omega*sin(Lb), v_n_eb(1)/(RN+hb);
           v_n_eb(2)*tan(Lb)/(RE+hb)+2*Earth_Omega*sin(Lb), (v_n_eb(1)*tan(Lb)+v_n_eb(3))/(RE+hb), v_n_eb(2)/(RE+hb)+2*Earth_Omega*cos(Lb);
           -2*v_n_eb(1)/(RN+hb), -2*v_n_eb(2)/(RE+hb)-2*Earth_Omega*cos(Lb), 0];
    F23 = [-(v_n_eb(2)^2)*(sec(Lb))^2/(RE+hb)-2*v_n_eb(2)*Earth_Omega*cos(Lb), 0, v_n_eb(2)^2*tan(Lb)/(RE+hb)^2-v_n_eb(1)*v_n_eb(3)/(RN+hb)^2;
           v_n_eb(1)*v_n_eb(2)*sec(Lb)^2/(RE+hb)+2*v_n_eb(1)*Earth_Omega*cos(Lb)-2*v_n_eb(3)*Earth_Omega*sin(Lb), 0, -(v_n_eb(1)*v_n_eb(2)*tan(Lb)+v_n_eb(2)*v_n_eb(3))/(RE+hb)^2;
           2*v_n_eb(2)*Earth_Omega*sin(Lb), 0, v_n_eb(2)^2/(RE+hb)^2+v_n_eb(1)^2/(RN+hb)^2-2*g0/r_e_eS];
    F32 = [1/(RN+hb), 0, 0;
           0, 1/((RE+hb)*cos(Lb)), 0;
           0, 0, -1];
    F33 = [0, 0, -v_n_eb(1)/(RN+hb)^2;
           v_n_eb(2)*sin(Lb)/((RE+hb)*cos(Lb)^2), 0, -v_n_eb(2)/(RE+hb)^2/cos(Lb);
           0, 0, 0];
    
    zr3 = zeros(3,3);
    F = [F11, F12, F13, zr3, C_n_b;
         F21, F22, F23, C_n_b, zr3;
         zr3, F32, F33, zr3, zr3;
         zr3, zr3, zr3, zr3, zr3;
         zr3, zr3, zr3, zr3, zr3];
    A = eye(15) + F*dt;
    SLy = 1e0;
    Sp = [SLy, 0, 0;
          0, SLy, 0;
          0, 0, 1];
    H = [zr3, zr3, -Sp, zr3, zr3;
         zr3, -eye(3), zr3, zr3, zr3];
     
    xp = A*x;
    Pp = A*PP*A' + Q;
    K = Pp*H'/(H*Pp*H' + R);
    x_k1 = xp + K*(z - H*xp);
    P_k1 = Pp - K*H*Pp;
    
end