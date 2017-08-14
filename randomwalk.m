clc;
close all;
clear;
% 
% v = VideoWriter('randomwalk.avi', 'Uncompressed AVI');
% open(v);

N = 100;
dt = 0.01;
t = 0:dt:100;
w = zeros(size(t));%10*sin(2*pi*t);
figure(1);

th = zeros(size(w));
    for i = 1:size(w,2)
       th(i+1) = th(i) + w(i)*dt;
    end
plot(t,th(1:end-1), 'LineWidth', 3.0, 'Color', 'blue');
hold on;

% frame = getframe(1);
% for i=1:10
%     writeVideo(v, frame.cdata);
% end
%pause(1);

for j=1:N
for num=1:N
    e = randn(size(w));
    w_hat = w + e;

    % subplot(1,2,1);
    % plot(t,w);
    % hold on;
    % plot(t,w_hat);
    % hold off;
    
    th_hat = zeros(size(w));
    for i = 1:size(w,2)
       th(i+1) = th(i) + w(i)*dt;
       th_hat(i+1) = th_hat(i) + w_hat(i)*dt;
    end

    % subplot(1,2,2);
    plot(t,th_hat(1:end-1));
    axis([0, N, -5, 5]);
    hold on;
    %hold off;
%     frame = getframe(1);
%     writeVideo(v, frame.cdata);
    %
end
pause(0.1)
end

plot(t,th(1:end-1), 'LineWidth', 3.0, 'Color', 'blue');
hold on;
plot(t,sqrt(t), 'LineWidth', 3.0, 'Color', 'red');
hold on;
% frame = getframe(1);
% writeVideo(v, frame.cdata);
% close(v);
