function result = IntKal(z)
    persistent A H Q R x P firstRun
    if isempty(firstRun)
        dt = 0.1;
        firstRun = 1;
        x = [0, 0]';
        P = 5*eye(2);
        A = [1 dt;
            0 1];
        H = [0, 1];
        Q = [1 0;
            0 3];
        R = 10;
    end
    
    xp = A*x;
    Pp = A*P*A' + Q;
    K = Pp*H'/(H*Pp*H' + R);
    x = xp + K*(z - H*xp);
    P = Pp - K*H*Pp;
    
    result = x;
end