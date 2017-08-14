function result = EulerGyro(prev, gy, dt)
    prevPhi = prev(1);
    prevTheta = prev(2);
    prevPsi = prev(3);
    p = gy(1); q = gy(2); r = gy(3);
    
    
    sinPhi = sin(prevPhi);
    cosPhi = cos(prevPhi);
    cosTheta = cos(prevTheta);
    tanTheta = tan(prevTheta);
    
    A = [1, sinPhi*tanTheta, cosPhi*tanTheta;
        0, cosPhi, -sinPhi;
        0, sinPhi/cosTheta, cosPhi/cosTheta];
    
    ang = [prevPhi, prevTheta, prevPsi]' + A*[p; q; r]*dt;
    
    result = ang;
end
