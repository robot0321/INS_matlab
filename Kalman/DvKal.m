function result = DvKal(z)

    persistent A H Q R x P firstRun

    if isempty(firstRun)
        firstRun=1;
        dt = 0.1;
        A = [1 dt;
            0 1];
        H = [1 0];
        Q = [1 0;
            0 3];
        R = 10;
        x = [0 0]';
        P = 5*eye(2);
    end

    xp = A*x;
    Pp = A*P*A' + Q;
    
    K = Pp*H'/(H*Pp*H' + R);
    x = xp + K*(z - H*xp);
    P = Pp - K*H*Pp;
    
    result = x;
end
