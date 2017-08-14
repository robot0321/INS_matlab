close all;
clear all;
clc;

figure();
hold on;
dt = 0.01;
N=1000;
t = dt:dt:N*dt;
for i=1:N
    a=zeros(1,N);
    for j=1:N-1
        a(j+1) = a(j) + 0.1*randn;
    end
    Y = fft(a);
    if (i==1)
        sum = abs(Y);
    else
        sum = (i-1)/i*sum + abs(Y)/i;
    end
    plot(t,a);
end
plot(t,sqrt(10*t),'LineWidth', 2,'Color','blue');
