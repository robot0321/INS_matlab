close all;
clear all;

dt=1;
t = 0:dt:100;
N = length(t);

pos = zeros(N,2);
vel = ones(N,2);

for i=1:N
    vel(i+1,:) = vel(i,:) + 1*[randn, randn];
    pos(i+1,:) = pos(i,:) + vel(i,:)*dt;    
end

%%

p = zeros(N,2);
for i=1:N
    v = randn(2,1);
    p(i,:) = TrackKalman(pos(i,:)'+v, dt)';
end

figure();
plot(p(:,1), p(:,2),'-o');
hold on;
plot(pos(:,1), pos(:,2),':*');