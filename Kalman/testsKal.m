clear all;
close all;

dt = 0.01;
t = 0:dt:10;

a = 3*sin(t) + 2*sin(3*t) + sin(10*t) + 10; 

N = length(t);
Xs = zeros(N,2);
Zs = zeros(N,1);


pos = 0; vel=80;

for k=1:N
    w = 1*randn;
    v = .5*randn;
    Zs(k) = a(k) + v;
    Xs(k,:) = DvKal(Zs(k));
  
end



figure();
plot(t, Xs(:,1), 'o-');
hold on;
plot(t, Zs, 'r.-');
plot(t, a, 'g:*');
hold off;

figure();
plot(t, Xs(:,2), 'o-');

% figure();
% plot(t, Xs(:,3), 'o-');
