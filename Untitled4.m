clear all;
close all;

M = importdata('hi.txt');

t = M(:,1);
acc = M(:,2:4);

ax = M(:,2);
ay = M(:,3);
az = M(:,4);
gx = M(:,5) - mean(M(1:490,5));
gy = M(:,6) - mean(M(1:490,6));
gz = M(:,7) - mean(M(1:490,7));
N = length(t);

mag = zeros(N,1);
for i=1:N
    mag(i) = norm(acc(i,:));
end

ax = ax./mag;
ay = ay./mag;
az = az./mag;


angA = zeros(N,2);
for i=1:N
    angA(i,:) = EulerAccel(ax(i),az(i),1);
end
Ca = angle2dcm(zeros(N,1), angA(:,1), angA(:,2));

qG = zeros(4,N);
angG = zeros(N,3);
qG(:,1) = angle2quat(0, angA(1,1), angA(1,2));
angG(1,:) = [angA(1,2), angA(1,1), 0]; % Y P R
for i=1:N-1
    qG(:,i+1) = QuatGyro(qG(:,i), [gx(i), gy(i), gz(i)], (M(i+1,1) - M(i,1))/1000);
    angG(i+1,:) = EulerGyro(angG(i,:), [gx(i), gy(i), gz(i)], (M(i+1,1) - M(i,1))/1000);
end
Cg = quat2dcm(qG');
Cgg = angle2dcm(angG(:,3), angG(:,2), angG(:,1));

attitude_visualize(Cg,Cgg);
