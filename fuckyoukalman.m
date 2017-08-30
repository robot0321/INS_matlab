quat_inverseINS;
clc;
clear;

load Qbody_data.mat;

dt = 0.01;
Earth_Omega = 7.292115e-5;
Earth_R_long = 6378137.0;
e = 0.0818191908425;
    
x_k = 0*ones(15,1);
P_k = 100*eye(15,15);

%%
P = init_P' - [0;0;Earth_R_long

q = init_q;

