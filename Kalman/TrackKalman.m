function result = TrackKalman(z,dt)
    persistent A H Q R x P firstRun
    
    if isempty(firstRun)
        A = eye(4,4); A(1,2) = dt; A(3,4) = dt;
        H = [1, 0, 0, 0;
            0, 0, 1, 0];
        Q = 1.0*eye(4);
        R = 50*eye(2);
        x = zeros(4,1);
        P = 100*eye(4);
        firstRun = 1;
    end
    
    xp = A*x;
    Pp = A*P*A' + Q;
    K = Pp*H'/(H*Pp*H' + R);
    x = xp + K*(z - H*xp);
    P = Pp - K*H*Pp;
    
    result = [x(1), x(3)];   
end