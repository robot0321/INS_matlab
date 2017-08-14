function result=QuatGyro(prev, gy, dt)
    A = [0,   -gy(1), -gy(2), -gy(3);
         gy(1),  0,   gy(3), -gy(2);
         gy(2), -gy(3),  0,  gy(1);
         gy(3), gy(2), -gy(1),  0];
    dq = 0.5*A*prev;
    result = prev + dq*dt;
end