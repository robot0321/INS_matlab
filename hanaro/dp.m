%data processing
clc;
close all;
clear;

M = importdata('realdata.txt');
MR = M(1,:);
for i=1:size(M,1)
    if(M(i,1)~=MR(end,1))
        MR(end+1,:) = M(i,:);
    end
end

% %sampling rate
% for i=1:size(MR,1)-1
%     a(i) = MR(i+1,1) - MR(i,1);
% end

%% altitude

alt_max = -100;
for i=1:size(MR,1)
   if(MR(i,16)~=1)
       rocket_start = i;
   end
   if(MR(i,14)==0)
       phase01 = i;
   end
   if(MR(i,14)~=2)
       phase12 = i;
   end
   if(MR(i,11) > alt_max)
      alt_max = MR(i,11); 
      alt_max_idx = i;
   end
end
alt_avg = mean(MR(1:490,11));

figure(1);
plot(MR(1:rocket_start,1)/1000, MR(1:rocket_start,11) - alt_avg);
hold on;
plot(MR(rocket_start:phase01,1)/1000, MR(rocket_start:phase01,11) - alt_avg);
hold on;
plot(MR(phase01:phase12,1)/1000, MR(phase01:phase12,11) - alt_avg);
hold on;
plot(MR(phase12:end,1)/1000, MR(phase12:end,11) - alt_avg);
hold off;
axis([0,85,-20,300]);
legend('ground', 'phase0', 'phase1', 'phase2','Location','northwest');
text(MR(rocket_start,1)/1000, MR(rocket_start,11) - alt_avg,'launch', 'horizontalAlignment', 'center', 'verticalAlignment', 'top');
text(MR(alt_max_idx,1)/1000, MR(alt_max_idx,11) - alt_avg,strcat('peak(',num2str(MR(alt_max_idx,11)- alt_avg),'m)'), 'horizontalAlignment', 'center', 'verticalAlignment', 'bottom');

text(MR(phase01,1)/1000, MR(phase01,11) - alt_avg,strcat('drogue(',num2str(MR(phase01,11)- alt_avg),'m, ', num2str(MR(phase01,11)-MR(alt_max_idx,11)),'m from peak)'), 'horizontalAlignment', 'left', 'verticalAlignment', 'middle');
text(MR(phase12,1)/1000, MR(phase12,11) - alt_avg,strcat('main(',num2str(MR(phase12,11)- alt_avg),'m, ', num2str(MR(phase12,11)-MR(alt_max_idx,11)),'m from peak)'), 'horizontalAlignment', 'left', 'verticalAlignment', 'middle');

N=1;
MR_avgN = zeros(size(MR,1)-N,1);
for i=1:size(MR,1)-N
    MR_avgN(i) = mean(MR(i:i+N-1,11));
end
%hold on;
%plot(MR(1:end-N,1),MR_avgN);

figure(2);

apx_vel = zeros(size(MR_avgN,1)-1,1);
for i=1:size(MR_avgN)-1
   apx_vel(i) = (MR_avgN(i+1) - MR_avgN(i))/(MR(i+1,1) - MR(i,1))*1000; 
end
plot(MR(1:end-N-1,1)/1000, apx_vel);

M=1;
vel_avgN = zeros(size(apx_vel,1)-M,1);
for i=1:size(apx_vel,1)-M
    vel_avgN(i) = mean(apx_vel(i:i+M-1));
end
apx_acc = zeros(size(vel_avgN,1)-1,1);
for i=1:size(vel_avgN)-1
   apx_acc(i) = (vel_avgN(i+1) - vel_avgN(i))/(MR(i+1,1) - MR(i,1))*1000; 
end

figure(3);
plot(MR(1:end-M-N-2,1)/1000, apx_acc/9.8);

%% attitude
ax = MR(:,2); ay = MR(:,3); az = MR(:,4);
gx = MR(:,5); gy = MR(:,6); gz = MR(:,7);
mx = MR(:,8); my = MR(:,9); mz = MR(:,10);

f_b = MR(:,2:4)';
w_b_ib = [gx'-mean(gx(1:490)); gy'-mean(gy(1:490)); gz'-mean(gz(1:490))];
mag = MR(:,8:10);

% load f_b, w_b_ib, init_P, init_C_n_b, init_V_n

dt = 0.1; %100Hz
Earth_Omega = 7.292115e-5;
Earth_R_short = 6356752.3142;
Earth_R_long = 6378137.0;
[ex, ey, ez] = ellipsoid(0,0,0,Earth_R_long,Earth_R_long,Earth_R_short,20);

% 초기값
P = [34.608270*pi/180, 127.204635*pi/180, R_surface(Earth_R_long, Earth_R_short, 34.608270*pi/180)]'; % 위도(L), 경도(l), 중심으로부터 거리(r)
avgf = [mean(f_b(1,1:490)),mean(f_b(2,1:490)),mean(f_b(3,1:490))];
init_rpy = acc2att(avgf);

qC = angle2dcm(init_rpy(1), init_rpy(2), init_rpy(3));
q = dcm2quat(qC)';
V_n = [0, 0, 0]';

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

figure(4);
subplot(3,1,1);
hold on;
plot(P(1,1:end-1)*180/pi)
subplot(3,1,2);
hold on;
plot(P(2,1:end-1)*180/pi)
subplot(3,1,3);
hold on;
plot(P(3,1:end-1))

figure(5);
%surf(ex, ey, ez);
%axis equal
hold on;
plot3(P(3,1:end-1).*cos(P(1,1:end-1)).*cos(P(2,1:end-1)), P(3,1:end-1).*cos(P(1,1:end-1)).*sin(P(2,1:end-1)), P(3,1:end-1).*sin(P(1,1:end-1)), 'LineWidth', 2, 'Marker', '.');
figure(6);
hold on;
plot(P(1,1:end-1)*180/pi,P(2,1:end-1)*180/pi);

figure();
attitude_visualize(qC,qC);
