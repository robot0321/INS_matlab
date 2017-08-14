function result = EulerKalman(gy,q,dt)
    
    persistent H Q R x P firstRun
    z = q;
    
    if isempty(firstRun)
       H = eye(4);
       Q = 1*eye(4);
       R = 1*eye(4);
       x = [1; 0; 0; 0];
       P = 1*eye(4);
       
       firstRun=1;
    end

    A = eye(4) + 0.5*dt*[0,   -gy(1), -gy(2), -gy(3);
                         gy(1),  0,   gy(3), -gy(2);
                         gy(2), -gy(3),  0,  gy(1);
                         gy(3), gy(2), -gy(1),  0];
    
    xp = A*x;
    Pp = A*P*A' + Q;
    K = Pp*H'/(H*Pp*H' + R);
    x = xp + K*(z - H*xp);
    P = Pp - K*H*Pp;
    
    result = x;
end